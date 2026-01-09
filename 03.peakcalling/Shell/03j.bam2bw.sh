#!/bin/bash

# copy from Jie
### genearte bw with RPKM norm. Add more options later

usage() { 
cat <<EOF
Usage: 12.bamCoverage.sh [-h] [-i input.bam] [-o out_dir] [-n name] [-s smooth]
Description: bam to bigwig. --binsize 50 -e 200 (--smoothLength 5*binsize)
Options:    
    -h           Print help and exit
    -i           input bam 
    -t           sort bam
    -n           bigwig initial
    -o           output path
    -s           smooth (--smoothLength 5*binsize) or not? T/F
EOF
    exit 1
}

smooth="F"
sort="T"
lens=250

while getopts ":i:t:o:n:s:l:h" flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        t) sort=${OPTARG};;
        o) out_path=${OPTARG};;
		n) name=${OPTARG};;
        s) smooth=${OPTARG};;
        l) lens=${OPTARG};;
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

if [[ ${smooth} = "T" ]]
then
    bamCoverage -b ${input} -o ${out_path}/${name}.bw --normalizeUsing RPKM --binSize 50 -e 200 --smoothLength ${lens} -p 4
else
    bamCoverage -b ${input} -o ${out_path}/${name}.bw --normalizeUsing RPKM --binSize 50 -e 200 -p 4
fi
#bamCoverage -b ${input} -o ${out_path}/${name}.bw --normalizeUsing RPKM # this is for RNA


