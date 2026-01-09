#!/bin/bash
### genearte bw from bam.

usage() { 
cat <<EOF
Usage: 02q.RNA.bam2bw.sh [-h] [-i input.bam] [-o out_dir] [-n name] [-m normalize]
Description: bam to bigwig.
Options:    
    -h           Print help and exit
    -i           input bam 
    -t           sort bam
    -n           bigwig initial
    -o           output path
    -m           normalize
EOF
    exit 1
}

sort="T"

while getopts ":i:to:n:m:h" flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        t) sort=${OPTARG};;
        o) out_path=${OPTARG};;
		n) name=${OPTARG};;
        m) normalize=${OPTARG};;
		h) usage;;
    esac
done

if [ -z "${input}" ] || [ -z "${out_path}" ] || [ -z "${name}" ]; then
    usage
fi

in_dir=`dirname $input`

if [[ ${sort} = "T" ]]; then
samtools sort -@ 2 -T ${out_path} ${input} > ${in_dir}/${name}.sorted.bam
if [[ -f ${in_dir}/${name}.sorted.bam ]] && [[ -s ${in_dir}/${name}.sorted.bam ]]; then
    mv ${in_dir}/${name}.sorted.bam ${input}
    echo "Input file has been sorted."
else
    echo "Failed to sort file."
    exit 1
fi
fi

samtools index ${input}

if [[ -z "${normalize}" ]]
then
    bamCoverage -b ${input} -o ${out_path}/${name}.bw -p 4
else
    bamCoverage -b ${input} -o ${out_path}/${name}.bw --normalizeUsing ${normalize} -p 4
fi
#bamCoverage -b ${input} -o ${out_path}/${name}.bw --normalizeUsing RPKM # this is for RNA


