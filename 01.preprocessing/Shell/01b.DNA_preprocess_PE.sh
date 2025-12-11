#!/usr/bin/bash
# example: bash 01b.DNA_preprocess_PE.sh -s XJ548 -g mm10 -m PT3 -p example-folder/

tscc=/tscc/projects
reachtool=${tscc}/package/Paired-Tag-master/reachtools/reachtools
rm_dup=rmdup_nh

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
    t)
      tscc="$OPTARG"
      echo "Processing option 't' with '${OPTARG}' argument"
      ;;
    l)
      reachtool="$OPTARG"
      echo "Processing option 'l' with '${OPTARG}' argument"
      ;;
    r)
      rm_dup="$OPTARG"
      echo "Processing option 'r' with '${OPTARG}' argument"
      ;;
    p)
      path="$OPTARG"
      echo "Processing option 'p' with '${OPTARG}' argument"
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

bash rename.sh ${sampleID}

fastqc_dir="${path}/00.fastqc"
fastq_dir="${path}/01.rawdata"
old_fastq_dir="${path}"
trim_dir="${path}/02.trimmed"
map_dir="${path}/03.mapping"
mtx_dir="${path}/04.matrices"
log_dir="${path}/log"

PT2="${tscc}/references/cell_id"
PT3="${tscc}/references/cell_id_407"
PT48="${tscc}/references/PairedTag48_384"

mm10_bowtie2="${tscc}/bowtie2_indexes/mm10"
hg38_bowtie2="${tscc}/bowtie2_indexes/hg38"
mix_bowtie2="${tscc}/genome_ref/GRCh38_mm10"

mm10_5k="${tscc}/reference/bin_ref/mm10.bin5k.txt"
hg38_5k="${tscc}/genome_ref/Paired-Tag/hg38/hg38.bin5k.txt"
mix_5k="${tscc}/genome_ref/Paired-Tag/mix/mix.bin5k.txt"

if [ -z ${genome+x} ]; then echo "genome is not specified. Default to mm10"; genome="mm10"; fi
if [ $genome == "mm10" ]; then ref=${mm10_bowtie2}; bin=${mm10_5k}; fi
if [ $genome == "hg38" ]; then ref=${hg38_bowtie2}; bin=${hg38_5k}; fi
if [ $genome == "mix" ]; then ref=${mix_bowtie2}; bin=${mix_5k}; fi

echo "barcoding mode: "${mode}
if [ $mode == "PT2" ]; then PT=${PT2}; combine="combine2"; fi
if [ $mode == "PT3" ]; then PT=${PT3}; combine="combine3"; fi
if [ $mode == "PT48" ]; then PT=${PT48}; combine="combine48"; fi

echo "Step00: FastQC!"
fastqc -t 8 -o ${fastqc_dir} ${fastq_dir}/${sampleID}_R?.fq.gz

echo "Step01: fastq combine!"
### Split fastq (Output: ${fastq_dir}/${sampleID}_R2_split_R2.fq.gz; ${fastq_dir}/${sampleID}_R2_split_R3.fq.gz)
echo "split long read2 into dir ${fastq_dir}/${sampleID}_split ...", 
${reachtool} splitR2 ${fastq_dir}/${sampleID}_R2.fq.gz 

mkdir -p ${fastq_dir}/${sampleID}_split/R1 ${fastq_dir}/${sampleID}_split/R3

ln -sf ${fastq_dir}/${sampleID}_R1.fq.gz ${fastq_dir}/${sampleID}_split/R1/${sampleID}_R1.fq.gz
ln -sf ${fastq_dir}/${sampleID}_R2_split_R2.fq.gz ${fastq_dir}/${sampleID}_split/R1/${sampleID}_R2.fq.gz
(${reachtool} ${combine} ${fastq_dir}/${sampleID}_split/R1/${sampleID}) 2>&1> ${log_dir}/${sampleID}_R1_qc.log 
ln -sf ${fastq_dir}/${sampleID}_split/R1/${sampleID}_DNA_combined.fq.gz ${fastq_dir}/${sampleID}_combined_R1.fq.gz

ln -sf ${fastq_dir}/${sampleID}_R2_split_R3.fq.gz ${fastq_dir}/${sampleID}_split/R3/${sampleID}_R1.fq.gz
ln -sf ${fastq_dir}/${sampleID}_R2_split_R2.fq.gz ${fastq_dir}/${sampleID}_split/R3/${sampleID}_R2.fq.gz
(${reachtool} ${combine} ${fastq_dir}/${sampleID}_split/R3/${sampleID}) 2>&1> ${log_dir}/${sampleID}_R3_qc.log 
ln -sf ${fastq_dir}/${sampleID}_split/R3/${sampleID}_DNA_combined.fq.gz ${fastq_dir}/${sampleID}_combined_R3.fq.gz

zcat ${fastq_dir}/${sampleID}_combined_R1.fq.gz | bowtie ${PT} - --norc -m 1 -v 1 -S ${fastq_dir}/${sampleID}_R1_BC.sam
zcat ${fastq_dir}/${sampleID}_combined_R3.fq.gz | bowtie ${PT} - --norc -m 1 -v 1 -S ${fastq_dir}/${sampleID}_R3_BC.sam
(${reachtool} convert2 ${fastq_dir}/${sampleID}_R1_BC.sam) 2>&1>> ${log_dir}/${sampleID}_R1_qc.log
(${reachtool} convert2 ${fastq_dir}/${sampleID}_R3_BC.sam) 2>&1>> ${log_dir}/${sampleID}_R3_qc.log

ln -sf ${log_dir}/${sampleID}_R1_qc.log ${log_dir}/${sampleID}_qc.log

if [[ -f "${fastq_dir}/${sampleID}_R1_BC_cov.fq.gz" && -f "${fastq_dir}/${sampleID}_R3_BC_cov.fq.gz" ]]
then 
	echo ${sampleID}" has been processed."
	rm ${fastq_dir}/${sampleID}_R1_BC.sam ${fastq_dir}/${sampleID}_R3_BC.sam
fi

# whether do merge or not?
if [[ -f "${old_fastq_dir}/${sampleID}_R1_BC_cov.fq.gz" ]]
then
    echo "previous sequencing file in "${old_fastq_dir}" has been found. Merge files for processing..."
    cat ${fastq_dir}/${sampleID}_R1_BC_cov.fq.gz ${old_fastq_dir}/${sampleID}_R1_BC_cov.fq.gz > ${fastq_dir}/${sampleID}_merge_R1_BC_cov.fq.gz
    cat ${fastq_dir}/${sampleID}_R_BC_cov.fq.gz ${old_fastq_dir}/${sampleID}_R3_BC_cov.fq.gz > ${fastq_dir}/${sampleID}_merge_R3_BC_cov.fq.gz
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
trim_galore --paired ${fastq_dir}/${salll}_R1_BC_cov.fq.gz ${fastq_dir}/${salll}_R3_BC_cov.fq.gz -o ${trim_dir} > ${trim_dir}/${salll}_PE.log 2>&1

echo "Step03: mapping!"
### need to retain all reads
(bowtie2 -x ${ref} -1 ${trim_dir}/${salll}_R1_BC_cov_val_1.fq.gz -2 ${trim_dir}/${salll}_R3_BC_cov_val_2.fq.gz -p 16 -S ${map_dir}/${salll}_${genome}.sam) 2>${map_dir}/${salll}.log
samtools sort -@ 16 -T ${map_dir} -o ${map_dir}/${salll}_${genome}_sorted.bam ${map_dir}/${salll}_${genome}.sam

if [[ -f "${map_dir}/${salll}_${genome}_sorted.bam" ]]
then 
	rm ${map_dir}/${salll}_${genome}.sam
fi

### bam2mtx
echo "Step04: matrix!"
### rmdup2 doesn't work on Paired-end file.
echo "extract R1 using sam flag 64..."
samtools view -@ 16 -h -f 64 ${map_dir}/${salll}_${genome}_sorted.bam | samtools view -@ 16 -bS > ${map_dir}/${salll}_${genome}_sorted_R1.bam ### extract mapped reads + first in pairs

echo "remove duplicates using reachtools..."
${reachtool} ${rm_dup} ${map_dir}/${salll}_${genome}_sorted_R1.bam
samtools view ${map_dir}/${salll}_${genome}_sorted_R1_rmdup.bam | awk '{print $1}' > ${map_dir}/${salll}_${genome}_sorted_R1_rmdup.readname

echo "extract dedup reads from R3..."
samtools view -H ${map_dir}/${salll}_${genome}_sorted.bam > ${map_dir}/${salll}_${genome}.head
samtools view -h ${map_dir}/${salll}_${genome}_sorted.bam | grep -Ff ${map_dir}/${salll}_${genome}_sorted_R1_rmdup.readname > ${map_dir}/${salll}_${genome}.sam
cat ${map_dir}/${salll}_${genome}.head ${map_dir}/${salll}_${genome}.sam | samtools view -bS -o ${map_dir}/${salll}_${genome}_sorted_SErmdup.bam

if [[ -f "${map_dir}/${salll}_${genome}_sorted_SErmdup.bam" ]]
then 
  rm ${map_dir}/${salll}_${genome}_sorted_R1_rmdup.readname ${map_dir}/${salll}_${genome}_sorted_R1.bam ${map_dir}/${salll}_${genome}.head ${map_dir}/${salll}_${genome}.sam
  mv ${map_dir}/${salll}_${genome}_sorted_SErmdup.bam ${map_dir}/${salll}_${genome}_sorted_rmdup.bam
  echo "done!"
fi

${reachtool} bam2Mtx2 ${map_dir}/${salll}_${genome}_sorted_rmdup.bam ${bin}
cp -r ${map_dir}/${salll}_${genome}_sorted_rmdup_mtx2 ${mtx_dir}

echo "Job is done!"
