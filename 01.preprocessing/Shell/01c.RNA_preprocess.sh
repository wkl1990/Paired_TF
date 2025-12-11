#!/usr/bin/bash

# example: bash 01c.RNA_preprocess.sh -s XJ554 -g mm10 -m PT3 -p example-folder/

tscc=/tscc/projects
star=/tscc/softwares
reachtool=${tscc}/package/Paired-Tag-master/reachtools/reachtools

while getopts ':s:g:m:p:t:l:r:h' opt; do
  case "$opt" in
    s)
      sampleID="$OPTARG"
      echo "Processing option 's' with '${OPTARG}' argument"
      ;;
    g)
      genome="$OPTARG"
      echo "Processing option 'g' with '${OPTARG}' argument"
      ;;
    m)
      mode="$OPTARG"
      echo "Processing option 'm' with '${OPTARG}' argument"
      ;;      
    p)
      path="$OPTARG"
      echo "Processing option 'p' with '${OPTARG}' argument"
      ;;
    t)
      tscc="$OPTARG"
      echo "Processing option 't' with '${OPTARG}' argument"
      ;;
    l)
      reachtool="$OPTARG"
      echo "Processing option 'l' with '${OPTARG}' argument"
      ;;
    r)
      star="$OPTARG"
      echo "Processing option 'r' with '${OPTARG}' argument"
      ;;
    h)
      echo "Usage: $(basename $0) [-s sampleID] [-g genome] [-m mode] [-p path]"
      exit 0
      ;;
    :)
      echo -e "option requires an argument.\nUsage: $(basename $0) [-s sampleID] [-g genome] [-m mode] [-p path]"
      exit 1
      ;;
    ?)
      echo -e "Invalid command option.\nUsage: $(basename $0) [-s sampleID] [-g genome] [-m mode] [-p path]"
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"


if [ -d "$path" ]; then
  cd $path
else
    echo -e "Error: ${path} not exist!"
    exit 1
fi

if [ ! -d "00.fastqc" ]; then
  mkdir 00.fastqc
fi
if [ ! -d "02.trimmed" ]; then
  mkdir 02.trimmed
fi
if [ ! -d "03.mapping" ]; then
  mkdir 03.mapping
fi
if [ ! -d "04.matrices" ]; then
  mkdir 04.matrices
fi
if [ ! -d "log" ]; then
  mkdir log
fi
if [ ! -d "featurecounts" ]; then
  mkdir featurecounts
fi

bash rename.sh ${sampleID}

fastqc_dir="${path}/00.fastqc"
fastq_dir="${path}/01.rawdata"
old_fastq_dir="${path}"
trim_dir="${path}/02.trimmed"
map_dir="${path}/03.mapping"
mtx_dir="${path}/04.matrices"
log_dir="${path}/log"
featurecounts_dir="${path}/featurecounts"

PT2="${tscc}/references/cell_id"
PT3="${tscc}/references/cell_id_407"
PT48="${tscc}/references/PairedTag48_384"

mm10_rna="${tscc}/reference/gene_ref/mm10.PairedTag.txt" #position, ensembl name, gene name
hg38_rna="${tscc}/genome_ref/Paired-Tag/hg38/hg38.gcode.v38.txt"
mix_rna="${tscc}/genome_ref/Paired-Tag/mix/mix.Paired-Tag.txt"

mm10_STAR="${tscc}/reference/star"
hg38_STAR="${tscc}/genome_ref/hg38_STAR"
mix_STAR="${tscc}/genome_ref/GRCh38_and_mm10/star"

mm10_gtf="${tscc}/reference/gtf_ref/gencode.vM25.annotation.gtf"
hg38_gtf="${tscc}/genome_ref/gencode.vH35.annotation.gtf"
mix_gtf="${tscc}/genome_ref/GRCh38_and_mm10/genes/genes.gtf.gz"

if [ $genome == "mm10" ]; then ref=${mm10_STAR}; gene=${mm10_rna}; gtf=${mm10_gtf}; STAR="${star}/miniconda3/envs/PairedTag/bin/STAR"; fi 
if [ $genome == "hg38" ]; then ref=${hg38_STAR}; gene=${hg38_rna}; gtf=${hg38_gtf}; STAR="${star}/miniconda3/bin/STAR"; fi #2.7.4a
if [ $genome == "mix" ]; then ref=${mix_STAR}; gene=${mix_rna}; gtf=${mix_gtf}; STAR="${star}/miniconda3/envs/PairedTag/bin/STAR"; fi ### 10X reference is built with 2.7.1

echo "barcoding mode: "${mode}
if [ $mode == "PT2" ]; then PT=${PT2}; combine="combine2"; fi
if [ $mode == "PT3" ]; then PT=${PT3}; combine="combine3"; fi 

echo "Step01: fastq combine!" 
(${reachtool} ${combine} ${fastq_dir}/${sampleID}) 2>&1> ${log_dir}/${sampleID}_qc.log 
ln -sf ${fastq_dir}/${sampleID}_RNA_combined.fq.gz ${fastq_dir}/${sampleID}_combined.fq.gz
zcat ${fastq_dir}/${sampleID}_combined.fq.gz | bowtie ${PT} - --norc -m 1 -v 1 -S ${fastq_dir}/${sampleID}_BC.sam
(${reachtool} convert2 ${fastq_dir}/${sampleID}_BC.sam) 2>&1>> ${log_dir}/${sampleID}_qc.log

if [[ -f "${fastq_dir}/${sampleID}_BC_cov.fq.gz" ]]
then 
     echo ${sampleID}" has been processed."
     rm ${fastq_dir}/${sampleID}_BC.sam
fi

# whether do merge or not?
if [[ -f "${old_fastq_dir}/${sampleID}_BC_cov.fq.gz" ]]
then
     echo "previous sequencing file in "${old_fastq_dir}" has been found. Merge files for processing..."
     cat ${fastq_dir}/${sampleID}_BC_cov.fq.gz ${old_fastq_dir}/${sampleID}_BC_cov.fq.gz > ${fastq_dir}/${sampleID}_merge_BC_cov.fq.gz
     sall=${sampleID}_merge
else
     sall=${sampleID}
fi

# whether do subsample or not? default is no subsampling
if [ -z ${sample+x} ]
then
     echo "no fastq sub sample is done"
     salll=${sall}
else
     seqtk sample -s 123 ${fastq_dir}/${sall}_BC_cov.fq.gz ${sample} | gzip > ${fastq_dir}/${sall}_sample_BC_cov.fq.gz
     salll=${sall}_sample
fi

echo "Step02: trimming!"
trim_galore ${fastq_dir}/${salll}_BC_cov.fq.gz -o ${trim_dir}
trim_galore -a AAAAAAAAAAAAAAAACCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${trim_dir}/${salll}_BC_cov_trimmed.fq.gz -o ${trim_dir} ### trim oligo-dT primer
trim_galore -a CCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${trim_dir}/${salll}_BC_cov_trimmed_trimmed.fq.gz -o ${trim_dir} ## trim N6 primer

echo "Step03: mapping!"
${STAR} --runThreadN 8 --genomeDir ${ref} --readFilesIn ${trim_dir}/${salll}_BC_cov_trimmed_trimmed_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix ${map_dir}/${salll}_${genome}_ --outSAMtype BAM Unsorted --quantMode GeneCounts
samtools view -h -F 256 ${map_dir}/${salll}_${genome}_Aligned.out.bam -b > ${map_dir}/${salll}_clean.bam
samtools sort -@ 16 -T ${map_dir} -o ${map_dir}/${salll}_${genome}_sorted.bam ${map_dir}/${salll}_clean.bam

if [[ -f "${map_dir}/${salll}_${genome}_sorted.bam" ]] 
then 
     rm ${map_dir}/${salll}_clean.bam 
fi

echo "Step04: matrix!"
${reachtool} rmdup2 ${map_dir}/${salll}_${genome}_sorted.bam
${reachtool} bam2Mtx2 ${map_dir}/${salll}_${genome}_sorted_rmdup.bam ${gene}
cp -r ${map_dir}/${salll}_${genome}_sorted_rmdup_mtx2 ${mtx_dir}
featureCounts -O -T 1 -s 1 -a ${gtf} -o ${featurecounts_dir}/${salll}.counts.txt ${map_dir}/${salll}_${genome}_sorted_rmdup.bam

echo "Job is done!"