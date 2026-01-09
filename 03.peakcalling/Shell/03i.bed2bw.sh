#!/usr/bin/bash

genome=mm
chrom=/tscc/Reference/genome/mm10/mm10.chrom.sizes

while getopts ':f:s:g:m:b:c:h' opt; do
  case "$opt" in
    f)
      fragment="$OPTARG"
      echo "Processing option 'f' with '${OPTARG}' argument"
      ;;
    s)
      sampleID="$OPTARG"
      echo "Processing option 's' with '${OPTARG}' argument"
      ;;
    g)
      genome="$OPTARG"
      echo "Processing option 'g' with '${OPTARG}' argument"
      ;;
    m)
      macs_path="$OPTARG"
      echo "Processing option 'm' with '${OPTARG}' argument"
      ;;      
    b)
      bigwig_path="$OPTARG"
      echo "Processing option 'b' with '${OPTARG}' argument"
      ;;      
    c)
      chrom="$OPTARG"
      echo "Processing option 'c' with '${OPTARG}' argument"
      ;;      
    h)
      echo "Usage: $(basename $0) [-s sampleID] [-g genome] [-m mode] [-p path]"
      exit 0
      ;;
    :)
      echo -e "option requires an argument.\nUsage: $(basename $0) [-f fragment] [-s sampleID] [-g genome] [-m macs_path] [-b bigwig_path] [-c chrom]"
      exit 1
      ;;
    ?)
      echo -e "Invalid command option.\nUsage: $(basename $0) [-f fragment] [-s sampleID] [-g genome] [-m macs_path] [-b bigwig_path] [-c chrom]"
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

if [ ! -d ${macs_path} ]; then
  mkdir -p ${macs_path}
fi
if [ ! -d ${bigwig_path} ]; then
  mkdir -p ${bigwig_path}
fi

bedgraph=${macs_path}/${sampleID}_treat_pileup.bdg
sort_bedgraph=${macs_path}/${sampleID}_treat_pileup.srt.bdg
bigwig=${bigwig_path}/${sampleID}_treat_pileup.srt.bw
clip=${macs_path}/${sampleID}_treat_pileup.srt.bdg.clip
sort_clip=${macs_path}/${sampleID}_treat_pileup.srt.bdg.sort.clip
debug_clip=${macs_path}/${sampleID}_treat_pileup.srt.bdg.sort.clip.debug


macs2 callpeak -t ${fragment} -f BED -n ${sampleID} -g ${genome} -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir ${macs_path}
grep 'chr' ${bedgraph} | egrep -v "chrUn|random" | sort -k1,1 -k2,2n > ${sort_bedgraph}
/home/kaw033/bin/bedGraphToBigWig ${sort_bedgraph} ${chrom} ${bigwig}

# bigwig bug for extension reads (https://gist.github.com/taoliu/2469050)
# bigwig bug for overlap regions (https://groups.google.com/a/soe.ucsc.edu/g/genome/c/_UqYFuzBwL4)

if [[ -f ${bigwig} ]]; then
	echo "bigwig is exist."
	filesize=`ls -l ${bigwig} | awk '{ print $5 }'`
	if [ $filesize -lt 10240 ]; then
		echo "bigwig seems very small, regenerate it!"
	if [[ -f ${sort_bedgraph} ]]; then
		echo "bdg exist, will generate the bigwig"
		bedtools slop -i ${sort_bedgraph} -g ${chrom} -b 0 | bedClip stdin ${chrom} ${clip}
		LC_COLLATE=C sort -k1,1 -k2,2n ${clip} > ${sort_clip}
		cat ${sort_clip} | awk 'BEGIN{OFS="\t"}{if (NR>1 && prev_chr==$1 && prev_chr_e<=$2) {print $0}; prev_chr=$1; prev_chr_e=$3;}' > ${debug_clip}
		bedGraphToBigWig ${debug_clip} ${chrom} ${bigwig}
		rm -f ${clip} ${sort_clip} ${debug_clip}
	elif [[ ! -f ${sort_bedgraph} ]]; then
		echo "Error: no bdg file, please check the macs2!"
	fi
	fi
elif [[ ! -f ${bigwig} ]] && [[ -f ${sort_bedgraph} ]]; then
	echo "bdg exist but bigwig is not exist, will generate the bigwig"
	bedtools slop -i ${sort_bedgraph} -g ${chrom} -b 0 | bedClip stdin ${chrom} ${clip}
	LC_COLLATE=C sort -k1,1 -k2,2n ${clip} > ${sort_clip}
	cat ${sort_clip} | awk 'BEGIN{OFS="\t"}{if (NR>1 && prev_chr==$1 && prev_chr_e<=$2) {print $0}; prev_chr=$1; prev_chr_e=$3;}' > ${debug_clip}
	bedGraphToBigWig ${debug_clip} ${chrom} ${bigwig}
	rm -f ${clip} ${sort_clip} ${debug_clip}
elif [[ ! -f ${bigwig} ]] && [[ ! -f ${sort_bedgraph} ]]; then
	echo "Error: no bdg file, please check the macs2!"
fi

echo 'Job is done!'


