#bedtools bamtobed -i ${1} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | tee ${2}.tagAlign | gzip -cn > ${2}.tagAlign.gz
bedtools bamtobed -i ${1} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -nc > ${2}.tagAlign.gz
#zcat ${2}.tagAlign.gz | awk 'BEGIN{OFS=FS="\t"}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -nc > ${2}.shifted.tagAlign.gz
zcat -f ${2}.tagAlign.gz | awk 'BEGIN {OFS="\t"} {if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} if ($2 >= $3) {if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1}} print $0}' | gzip -nc > ${2}.shifted.tagAlign.gz

# PE
#bedtools bamtobed -bedpe -mate1 -i ${1} | gzip -nc > ${2}.bedpe.gz
#zcat -f {} | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${2}.tagAlign.gz
#zcat -f ${2}.tagAlign.gz | awk 'BEGIN {OFS="\t"} {if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} if ($2 >= $3) {if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1}} print $0}' | gzip -nc > ${2}.shifted.tagAlign.gz

prefix=`basename ${2}`
macs2 callpeak -t ${2}.shifted.tagAlign.gz -f BED -n ${prefix} -g mm -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir ${3}
