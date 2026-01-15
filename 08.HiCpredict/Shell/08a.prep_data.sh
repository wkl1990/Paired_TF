path=/tscc/projects/PairedTF
# Prepare marker gene bed 
cat /tscc/Reference/genes/mm10/modified_gencode.vM23.gene.bed | grep -E "Slc17a7$|Foxp2$|Calb1$|Slc1a3$|Rorb$|Sulf1$|Tcerg1l$|Tshz2$|Elfn1$|Vip$|Gad1$|Sst$|Plp1$|Pdgfra$|Flt1$|Cped1$" > ${path}/hic_pred/marker_genes.bed
# generate CTCF bw
cat ${path}/CTCF/celltype_subclass.txt | while read id 
do
celltype=`echo $id | awk '{FS="\t"}{print $2}'`
celltype_name=`echo $id | awk '{FS="\t"}{print $1}'`
subclass_name=`echo $id | awk '{FS="\t"}{print $4}' | sed 's/\./_/g;s/-/_/g'`
CTCF_bed=${path}/fragment/CTCF.${celltype_name}.bed.gz
CTCF_bam=${path}/bed2bam/CTCF.${celltype_name}.bam
CTCF_sort_bam=${path}/bed2bam/CTCF.${celltype_name}.sort.bam
CTCF_bigwig=${path}/bam2bigwig/CTCF.${celltype_name}.bigwig
if [[ -f ${CTCF_bed} ]] && [[ ! -f ${CTCF_bigwig} ]]; then
echo ${celltype_name}
bedtools bedtobam -i ${CTCF_bed} -g /tscc/Reference/genome/mm10/mm10.chrom.sizes > ${CTCF_bam}
samtools sort -o ${CTCF_sort_bam} ${CTCF_bam}
samtools index ${CTCF_sort_bam}
bamCoverage --normalizeUsing RPKM --binSize 1 --bam ${CTCF_sort_bam} -o ${CTCF_bigwig} &
else
echo "check ${celltype_name} CTCF file!"
fi
done

# generate ATAC bw
cat ${path}/CTCF/celltype_subclass.txt | while read id 
do
celltype=`echo $id | awk '{FS="\t"}{print $2}'`
celltype_name=`echo $id | awk '{FS="\t"}{print $1}'`
subclass_name=`echo $id | awk '{FS="\t"}{print $4}' | sed 's/\./_/g;s/-/_/g'`
ATAC_bam=/tscc/projects/CEMBA2/${subclass_name}.srt.bam
ATAC_bigwig=/tscc/projects/CEMBA2/${subclass_name}.srt.bigwig
if [[ -f ${ATAC_bam} ]] && [[ ! -f ${ATAC_bigwig} ]]; then
echo ${subclass_name}
bamCoverage --normalizeUsing RPKM --binSize 1 --bam ${ATAC_bam} -o ${ATAC_bigwig} &
else
echo "check ${subclass_name} ATAC file!"
fi
done


# normalize CTCF bw
cat ${path}/CTCF/celltype_subclass.txt | while read id 
do
celltype=`echo $id | awk '{FS="\t"}{print $2}'`
celltype_name=`echo $id | awk '{FS="\t"}{print $1}'`
subclass_name=`echo $id | awk '{FS="\t"}{print $4}' | sed 's/\./_/g;s/-/_/g'`
CTCF_bigwig=${path}/bam2bigwig/CTCF.${celltype_name}.bigwig
ATAC_bigwig=/tscc/projects/CEMBA2/${subclass_name}.srt.bigwig
CTCF_normalized_bigwig=${path}/bam2bigwig/CTCF.${celltype_name}.normalized.bigwig
if [[ -f ${CTCF_bigwig} ]] && [[ -f ${ATAC_bigwig} ]] && [[ ! -f ${CTCF_normalized_bigwig} ]]; then
echo ${celltype_name}
bigwigCompare --binSize 1 -b1 ${CTCF_bigwig} -b2 ${ATAC_bigwig} -o ${CTCF_normalized_bigwig} &
else
echo "check ${celltype_name} CTCF normalized bigwig!"
fi
done

# generate input data dir
cat ${path}/CTCF/celltype_subclass.txt | while read id 
do
celltype=`echo $id | awk '{FS="\t"}{print $2}'`
celltype_name=`echo $id | awk '{FS="\t"}{print $1}'`
subclass_name=`echo $id | awk '{FS="\t"}{print $4}' | sed 's/\./_/g;s/-/_/g'`
ATAC_bigwig=/tscc/projects/CEMBA2/${subclass_name}.srt.bigwig
CTCF_normalized_bigwig=${path}/bam2bigwig/CTCF.${celltype_name}.normalized.bigwig
bam2bw_dir=${path}/data/mm10/${celltype}_bam2bw
hic_dir=${path}/data/mm10/${celltype}/hic_matrix
if [[ -f ${CTCF_normalized_bigwig} ]] && [[ -f ${ATAC_bigwig} ]] && [[ ! -d ${bam2bw_dir} ]]; then
echo ${celltype_name}
mkdir -p ${bam2bw_dir}/genomic_features
ln -s ${CTCF_normalized_bigwig} ${bam2bw_dir}/genomic_features/ctcf_log2fc.bw
ln -s ${ATAC_bigwig} ${bam2bw_dir}/genomic_features/atac.bw
ln -s ${hic_dir} ${bam2bw_dir}
else
echo "check ${celltype_name} CTCF dir!"
fi
done

# generate control bw
cat ${path}/CTCF/celltype_subclass.txt | while read id 
do
celltype=`echo $id | awk '{FS="\t"}{print $2}'`
celltype_name=`echo $id | awk '{FS="\t"}{print $1}'`
subclass_name=`echo $id | awk '{FS="\t"}{print $4}' | sed 's/\./_/g;s/-/_/g'`
ATAC_bigwig=/tscc/projects/CEMBA2/${subclass_name}.srt.bigwig
CTCF_normalized_bigwig=${path}/bam2bigwig/CTCF.${celltype_name}.normalized.bigwig
hic_dir=${path}/data/mm10/${celltype}/hic_matrix
CTCF_normalized_bedgraph=${path}/bam2bigwig/CTCF.${celltype_name}.normalized.bedgraph
CTCF_zero_bedgraph=${path}/bam2bigwig/CTCF.${celltype_name}.zero.bedgraph
CTCF_zero_bigwig=${path}/bam2bigwig/CTCF.${celltype_name}.zero.bw
zero_dir=${path}/data/mm10/${celltype}_zeroCTCF
if [[ -f ${CTCF_normalized_bigwig} ]] && [[ ! -f ${CTCF_zero_bigwig} ]]; then
echo ${celltype_name}
bigWigToBedGraph ${CTCF_normalized_bigwig} ${CTCF_normalized_bedgraph}
cat ${CTCF_normalized_bedgraph} | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,0}' | sort -k1,1 -k2,2n > ${CTCF_zero_bedgraph}
bedGraphToBigWig ${CTCF_zero_bedgraph} /tscc/Reference/genome/mm10/mm10.chrom.sizes ${CTCF_zero_bigwig}
else
echo "check ${celltype_name} CTCF zero bigwig!"
fi
if [[ -f ${CTCF_zero_bigwig} ]] && [[ -f ${ATAC_bigwig} ]] && [[ ! -d ${zero_dir} ]]; then
echo ${celltype_name}
mkdir -p ${zero_dir}/genomic_features
ln -s ${CTCF_zero_bigwig} ${zero_dir}/genomic_features/ctcf_log2fc.bw
ln -s ${ATAC_bigwig} ${zero_dir}/genomic_features/atac.bw
ln -s ${hic_dir} ${zero_dir}
else
echo "check ${celltype_name} zero dir!"
fi
done


