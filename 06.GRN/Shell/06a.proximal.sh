# defined proximal binding within promoter (TSS +/- 1kb)
TF=$1

path=/tscc/projects/PairedTF
peak_file=${path}/${TF}.peaks.union.bed
ref_file=/tscc/projects/ref_genome/mm10-refseq.TSS.1kb.bed

if [ -f $peak_file ]
then
bedtools intersect -wao -f 0.5 -a ${peak_file} -b ${ref_file} > ${path}/${TF}.peaks.union.tssUpDn1k.bed
awk 'BEGIN{FS=OFS="\t"}{if($4!="."){print $1,$2,$3,"proximal",$8}else{print $1,$2,$3,"distal","nan"}}' ${path}/${TF}.peaks.union.tssUpDn1k.bed | sort -k1,1 -k2,2n | uniq >${path}/${TF}.peaks.union.tssUpDn1k.proximal.distal.bed
awk 'BEGIN{FS=OFS="\t"}{print $1"_"$2"_"$3,$4,$5}' ${path}/${TF}.peaks.union.tssUpDn1k.proximal.distal.bed | sort -k1,1 > ${path}/${TF}.peaks.union.tssUpDn1k.proximal.distal.peaks
sleep 1
else
echo "$peak_file is not exist!"
fi
