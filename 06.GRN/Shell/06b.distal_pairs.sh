# identify gene peak pairs within 1000Kb

TF=$1

path=/tscc/projects/PairedTF
peak_file=${path}/${TF}.peaks.union.bed
ref_file=/tscc/projects/reference/gtf_ref/gene.tssUpDn1000k.bed

if [ -f $peak_file ]
then
echo $TF
bedtools intersect -wao -a ${peak_file} -b ${ref_file} > ${path}/${TF}.peaks.union.tssUpDn1000k.bed
awk 'BEGIN{FS=OFS="\t"}{if($4!="."){print $10,$1":"$2"-"$3}}' ${path}/${TF}.peaks.union.tssUpDn1000k.bed | sort -k2,2 | uniq > ${path}/${TF}.peaks.union.tssUpDn1000k.gene.pairs
sleep 1
else
echo "$peak_file is not exist!"
fi

