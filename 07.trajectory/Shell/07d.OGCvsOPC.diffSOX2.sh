export PERL5LIB=~/bin/lib/site_perl/
export PATH=$PATH:~/bin/bin

# SOX2 OGC and OPC call peak
for celltype in OPC OGC
do
echo $celltype
nohup bash /path/07.trajectory/Shell/07c.bam2peak.sh ${path}/bam/SOX2_${celltype}.bam ${path}/OGCvsOPC/peaks/SOX2_${celltype} ${path}/OGCvsOPC/peaks/ 2>&1 > ${path}/OGCvsOPC/peaks/SOX2_${celltype}.bam2peak.log &
sleep 1
#loadavg
done

# merge SOX2 OGC and OPC peaks 
Rscript=/tscc/softwares/miniconda3/envs/scRNA/bin/Rscript
${Rscript} ${path}/07.trajectory/R/07c.peak_merge.R -i ${path}/OGCvsOPC/SOX2.naiveSummitList.list \
        -g mm10 \
        --extend 500 \
        --blacklist /tscc/reference/blacklist/mm10-blacklist.v2.bed3 \
        --chromSize /tscc/Reference/genome/mm10/mm10.chrom.sizes \
        -d ${path}/OGCvsOPC/peaks/ -o SOX2.OGCvsOPC.merge

cat ${path}/OGCvsOPC/peaks/SOX2.OGCvsOPC.merge.filteredNfixed.union.peakSet | grep -P 'chr[0-9XY]+(?!_)' | sed '1d' | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$7,$11}' | sort -k1,1 -k2,2n | uniq > ${path}/OGCvsOPC/peaks/SOX2.OGCvsOPC.merge.union.bed

# get matrix
multiBamSummary BED-file --BED ${path}/OGCvsOPC/peaks/SOX2.OGCvsOPC.merge.union.bed --bamfiles ${path}/bam/SOX2_OPC.bam ${path}/bam/SOX2_OGC.bam -p 2 --smartLabels -o ${path}/OGCvsOPC/mtx/SOX2_peak_q01.npz --outRawCounts ${path}/OGCvsOPC/mtx/SOX2_peak_q01.tab --scalingFactors ${path}/OGCvsOPC/mtx/SOX2_peak_q01.txt &

# diffreps
gunzip ${path}/OGCvsOPC/peaks/SOX2_OGC.shifted.tagAlign.gz 
gunzip ${path}/OGCvsOPC/peaks/SOX2_OPC.shifted.tagAlign.gz

nohup ~/bin/bin/diffReps.pl -tr ${path}/OGCvsOPC/peaks/SOX2_OGC.shifted.tagAlign -co ${path}/OGCvsOPC/peaks/SOX2_OPC.shifted.tagAlign -ch /tscc/Reference/genome/mm10/mm10.chrom.sizes -re ${path}/OGCvsOPC/SOX2_OGC_OPC_diffreps.gt.txt -me gt 2>&1 > ${path}/OGCvsOPC/SOX2_OGC_OPC_diffreps.gt.log &

cat ${path}/OGCvsOPC/SOX2_OGC_OPC_diffreps.gt.txt | grep "Up" > ${path}/OGCvsOPC/SOX2_OGC_OPC_diffreps.gt.pos.bed
cat ${path}/OGCvsOPC/SOX2_OGC_OPC_diffreps.gt.txt | grep "Down" > ${path}/OGCvsOPC/SOX2_OGC_OPC_diffreps.gt.neg.bed

# annotate diff peaks to proximal and distal
for diff in pos neg
do  
diff_file=${path}/OGCvsOPC/SOX2_peakq01_cpm_log2FC_${diff}.bed
bedtools intersect -wao -f 0.5 -a ${diff_file} -b /tscc/annot/gencode.vM23.gene.tssUpDn1k.bed > ${path}/OGCvsOPC/SOX2_peakq01_cpm_log2FC_${diff}.tssUpDn1k.bed
awk 'BEGIN{FS=OFS="\t"}{if($6!="."){print $1,$2,$3,$5,"proximal",$12}else{print $1,$2,$3,$5,"distal","nan"}}' ${path}/OGCvsOPC/SOX2_peakq01_cpm_log2FC_${diff}.tssUpDn1k.bed | sort -k1,1 -k2,2n | uniq >${path}/OGCvsOPC/SOX2_peakq01_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.bed
awk 'BEGIN{FS=OFS="\t"}{print $1"_"$2"_"$3,$4,$5,$6}' ${path}/OGCvsOPC/SOX2_peakq01_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.bed | sort -k1,1 > ${path}/OGCvsOPC/SOX2_peakq01_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.peaks
cat ${path}/OGCvsOPC/SOX2_peakq01_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.peaks | sed '1i peak\tlogFC\ttype\tgene' | awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}$3=="proximal"{print $0}' > ${path}/OGCvsOPC/SOX2_peakq01_cpm_log2FC_${diff}.proximal.peaks
done

# generate SOX2 OGC and OPC bigwig
for celltype in OPC OGC
do
echo $celltype
bdg=${path}/OGCvsOPC/peaks/SOX2_${celltype}_treat_pileup.bdg
bdg_srt=${path}/OGCvsOPC/peaks/SOX2_${celltype}_treat_pileup.srt.bdg
bigwig=${path}/OGCvsOPC/bw/SOX2_${celltype}_treat_pileup.srt.bw
clip=${path}/OGCvsOPC/peaks/SOX2_${celltype}_treat_pileup.srt.bdg.clip
clip_srt=${path}/OGCvsOPC/peaks/SOX2_${celltype}_treat_pileup.srt.bdg.sort.clip
clip_debug=${path}/OGCvsOPC/peaks/SOX2_${celltype}_treat_pileup.srt.bdg.sort.clip.debug
grep 'chr' ${bdg} | egrep -v "chrUn|random" | sort -k1,1 -k2,2n > ${bdg_srt}
/home/kaw033/bin/bedGraphToBigWig ${bdg_srt} /tscc/Reference/genome/mm10/mm10.chrom.sizes ${bigwig}
if [[ -f ${bigwig} ]]; then
	echo "bigwig is exist."
	filesize=`ls -l ${bigwig} | awk '{ print $5 }'`
	if [ $filesize -lt 10240 ]; then
		echo "bigwig seems very small, regenerate it!"
	if [[ -f ${bdg_srt} ]]; then
		echo "bdg exist, will generate the bigwig"
		bedtools slop -i ${bdg_srt} -g /tscc/Reference/genome/mm10/mm10.chrom.sizes -b 0 | bedClip stdin /tscc/Reference/genome/mm10/mm10.chrom.sizes ${clip}
		LC_COLLATE=C sort -k1,1 -k2,2n ${clip} > ${clip_srt}
		cat ${clip_srt} | awk 'BEGIN{OFS="\t"}{if (NR>1 && prev_chr==$1 && prev_chr_e<=$2) {print $0}; prev_chr=$1; prev_chr_e=$3;}' > ${clip_debug}
		bedGraphToBigWig ${clip_debug} /tscc/Reference/genome/mm10/mm10.chrom.sizes ${bigwig}
		rm -f ${clip} ${clip_srt} ${clip_debug}
	elif [[ ! -f ${bdg_srt} ]]; then
		echo "Error: no bdg file, please check the macs2!"
	fi
	fi
elif [[ ! -f ${bigwig} ]] && [[ -f ${bdg_srt} ]]; then
	echo "bdg exist but bigwig is not exist, will generate the bigwig"
	bedtools slop -i ${bdg_srt} -g /tscc/Reference/genome/mm10/mm10.chrom.sizes -b 0 | bedClip stdin /tscc/Reference/genome/mm10/mm10.chrom.sizes ${clip}
	LC_COLLATE=C sort -k1,1 -k2,2n ${clip} > ${clip_srt}
	cat ${clip_srt} | awk 'BEGIN{OFS="\t"}{if (NR>1 && prev_chr==$1 && prev_chr_e<=$2) {print $0}; prev_chr=$1; prev_chr_e=$3;}' > ${clip_debug}
	bedGraphToBigWig ${clip_debug} /tscc/Reference/genome/mm10/mm10.chrom.sizes ${bigwig}
	rm -f ${clip} ${clip_srt} ${clip_debug}
elif [[ ! -f ${bigwig} ]] && [[ ! -f ${bdg_srt} ]]; then
	echo "Error: no bdg file, please check the macs2!"
fi
echo 'Done'
done


