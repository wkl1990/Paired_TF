export PERL5LIB=~/bin/lib/site_perl/
export PATH=$PATH:~/bin/bin

# generate CTCF OGC and OPC bam 
for celltype in OGC OPC
do
echo $celltype
samtools merge -@ 16 ${path}/OGCvsOPC/bam/CTCF_${celltype}.bam ${path}/OGCvsOPC/scbam/metacell_${celltype}.MC-CTCF-RNPII_DNA_*.bam
samtools sort -@ 16 -T ${path}/OGCvsOPC/bam/ ${path}/OGCvsOPC/bam/CTCF_${celltype}.bam > ${path}/OGCvsOPC/bam/CTCF_${celltype}.sorted.bam
if [[ -f ${path}/OGCvsOPC/bam/CTCF_${celltype}.sorted.bam ]] && [[ -s ${path}/OGCvsOPC/bam/CTCF_${celltype}.sorted.bam ]]; then
    mv -f ${path}/OGCvsOPC/bam/CTCF_${celltype}.sorted.bam ${path}/OGCvsOPC/bam/CTCF_${celltype}.bam
    echo "${celltype} has been sorted."
else
    echo "Failed to sort ${celltype}."
fi
samtools index ${path}/OGCvsOPC/bam/CTCF_${celltype}.bam
sleep 1
#loadavg
done

# CTCF OGC and OPC call peak
ls ${path}/OGCvsOPC/bam/CTCF*bam | while read id 
do
file=$(basename $id)
name=${file%.*}
echo $name
nohup bash ${path}/07.trajectory/Shell/07c.bam2peak.sh ${id} ${path}/OGCvsOPC/peaks/${name} ${path}/OGCvsOPC/peaks/ 2>&1 > ${path}/OGCvsOPC/peaks/${name}.bam2peak.log &
sleep 3
#loadavg
done

# merge CTCF OGC and OPC peak
Rscript=/tscc/softwares/miniconda3/envs/scRNA/bin/Rscript
${Rscript} ${path}/07.trajectory/R/07c.peak_merge.R -i ${path}/OGCvsOPC/CTCF.naiveSummitList.list \
        -g mm10 \
        --extend 500 \
        --blacklist /tscc/reference/blacklist/mm10-blacklist.v2.bed3 \
        --chromSize /tscc/Reference/genome/mm10/mm10.chrom.sizes \
        -d ${path}/OGCvsOPC/peaks/ -o CTCF.OGCvsOPC.merge

cat ${path}/OGCvsOPC/peaks/CTCF.OGCvsOPC.merge.filteredNfixed.union.peakSet | grep -P 'chr[0-9XY]+(?!_)' | sed '1d' | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$7,$11}' | sort -k1,1 -k2,2n | uniq > ${path}/OGCvsOPC/peaks/CTCF.OGCvsOPC.merge.union.bed

# get matrix
multiBamSummary BED-file --BED ${path}/OGCvsOPC/peaks/CTCF.OGCvsOPC.merge.union.bed --bamfiles ${path}/OGCvsOPC/bam/CTCF*bam -p 2 --smartLabels -o ${path}/OGCvsOPC/mtx/CTCF_peak_q01.npz --outRawCounts ${path}/OGCvsOPC/mtx/CTCF_peak_q01.tab --scalingFactors ${path}/OGCvsOPC/mtx/CTCF_peak_q01.txt &

# diffreps
gunzip ${path}/OGCvsOPC/peaks/CTCF_OGC.shifted.tagAlign.gz 
gunzip ${path}/OGCvsOPC/peaks/CTCF_OPC.shifted.tagAlign.gz
nohup ~/bin/bin/diffReps.pl -tr ${path}/OGCvsOPC/peaks/CTCF_OGC.shifted.tagAlign -co ${path}/OGCvsOPC/peaks/CTCF_OPC.shifted.tagAlign -ch /tscc/Reference/genome/mm10/mm10.chrom.sizes -re ${path}/OGCvsOPC/CTCF_OGC_OPC_diffreps.gt.txt -me gt 2>&1 > ${path}/OGCvsOPC/CTCF_OGC_OPC_diffreps.gt.log &
cat ${path}/OGCvsOPC/CTCF_OGC_OPC_diffreps.gt.txt | grep "Up" > ${path}/OGCvsOPC/CTCF_OGC_OPC_diffreps.gt.pos.bed
cat ${path}/OGCvsOPC/CTCF_OGC_OPC_diffreps.gt.txt | grep "Down" > ${path}/OGCvsOPC/CTCF_OGC_OPC_diffreps.gt.neg.bed

# annotate diff peaks to proximal and distal
for diff in pos neg
do  
diff_file=${path}/OGCvsOPC/CTCF_peakq01_cpm_log2FC_${diff}.bed
bedtools intersect -wao -f 0.5 -a ${diff_file} -b /tscc/annot/gencode.vM23.gene.tssUpDn1k.bed >${path}/OGCvsOPC/CTCF_peakq01_cpm_log2FC_${diff}.tssUpDn1k.bed
awk 'BEGIN{FS=OFS="\t"}{if($6!="."){print $1,$2,$3,$5,"proximal",$12}else{print $1,$2,$3,$5,"distal","nan"}}' ${path}/OGCvsOPC/CTCF_peakq01_cpm_log2FC_${diff}.tssUpDn1k.bed | sort -k1,1 -k2,2n | uniq > ${path}/OGCvsOPC/CTCF_peakq01_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.bed
awk 'BEGIN{FS=OFS="\t"}{print $1"_"$2"_"$3,$4,$5,$6}' ${path}/OGCvsOPC/CTCF_peakq01_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.bed | sort -k1,1 > ${path}/OGCvsOPC/CTCF_peakq01_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.peaks
cat ${path}/OGCvsOPC/CTCF_peakq01_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.peaks | sed '1i peak\tlogFC\ttype\tgene' | awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}$3=="proximal"{print $0}' > ${path}/OGCvsOPC/CTCF_peakq01_cpm_log2FC_${diff}.proximal.peaks
done


# generate CTCF OGC and OPC bigwig
ls ${path}/OGCvsOPC/bam/CTCF*bam | while read id
do
file=$(basename $id) 
name=${file%.*}
echo $name
nohup bash ${path}/03.peakcalling/Shell/03j.bam2bw.sh -i ${id} -o ${path}/OGCvsOPC/bigwig/ -n ${name} -s T -t F 2>&1 > ${path}/OGCvsOPC/bigwig/${name}.bam2bw.log &
sleep 3
#loadavg
done

ls ${path}/OGCvsOPC/peaks/CTCF*_treat_pileup.bdg | while read id
do
file=$(basename $id)
celltype=`echo $file | sed 's/CTCF_//g;s/_treat_pileup.bdg//g'`
echo $celltype
grep 'chr' ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.bdg | egrep -v "chrUn|random" | sort -k1,1 -k2,2n > ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg
bedGraphToBigWig ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg /tscc/Reference/genome/mm10/mm10.chrom.sizes ${path}/OGCvsOPC/bigwig/CTCF_${celltype}_treat_pileup.srt.bw
if [[ -f ${path}/OGCvsOPC/bigwig/CTCF_${celltype}_treat_pileup.srt.bw ]]; then
    echo "bigwig is exist."
    filesize=`ls -l ${path}/OGCvsOPC/bigwig/CTCF_${celltype}_treat_pileup.srt.bw | awk '{ print $5 }'`
    if [ $filesize -lt 10240 ]; then
        echo "bigwig seems very small, regenerate it!"
    if [[ -f ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg ]]; then
        echo "bdg exist, will generate the bigwig"
        bedtools slop -i ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg -g /tscc/Reference/genome/mm10/mm10.chrom.sizes -b 0 | bedClip stdin /tscc/Reference/genome/mm10/mm10.chrom.sizes ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.clip
        LC_COLLATE=C sort -k1,1 -k2,2n ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.clip > ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip
        cat ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip | awk 'BEGIN{OFS="\t"}{if (NR>1 && prev_chr==$1 && prev_chr_e<=$2) {print $0}; prev_chr=$1; prev_chr_e=$3;}' > ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip.debug
        bedGraphToBigWig ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip.debug /tscc/Reference/genome/mm10/mm10.chrom.sizes ${path}/OGCvsOPC/bigwig/CTCF_${celltype}_treat_pileup.srt.bw
        rm -f ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.clip ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip.debug
    elif [[ ! -f ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg ]]; then
        echo "Error: no bdg file, please check the macs2!"
    fi
    fi
elif [[ ! -f ${path}/OGCvsOPC/bigwig/CTCF_${celltype}_treat_pileup.srt.bw ]] && [[ -f ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg ]]; then
    echo "bdg exist but bigwig is not exist, will generate the bigwig"
    bedtools slop -i ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg -g /tscc/Reference/genome/mm10/mm10.chrom.sizes -b 0 | bedClip stdin /tscc/Reference/genome/mm10/mm10.chrom.sizes ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.clip
    LC_COLLATE=C sort -k1,1 -k2,2n ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.clip > ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip
    cat ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip | awk 'BEGIN{OFS="\t"}{if (NR>1 && prev_chr==$1 && prev_chr_e<=$2) {print $0}; prev_chr=$1; prev_chr_e=$3;}' > ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip.debug
    bedGraphToBigWig ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip.debug /tscc/Reference/genome/mm10/mm10.chrom.sizes ${path}/OGCvsOPC/bigwig/CTCF_${celltype}_treat_pileup.srt.bw
    rm -f ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.clip ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg.sort.clip.debug
elif [[ ! -f ${path}/OGCvsOPC/bigwig/CTCF_${celltype}_treat_pileup.srt.bw ]] && [[ ! -f ${path}/OGCvsOPC/peaks/CTCF_${celltype}_treat_pileup.srt.bdg ]]; then
    echo "Error: no bdg file, please check the macs2!"
fi
echo 'Done'
done

