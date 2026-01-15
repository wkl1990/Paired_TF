path=/tscc/projects/PairedTF
# prediction marker regions
cat ${path}/CTCF/celltype_subclass.txt | while read id 
do
celltype=`echo $id | awk '{FS="\t"}{print $2}'`
celltype_name=`echo $id | awk '{FS="\t"}{print $1}'`
subclass_name=`echo $id | awk '{FS="\t"}{print $4}' | sed 's/\./_/g;s/-/_/g'`
bam2bw_dir=${path}/data/mm10/${celltype}_bam2bw
bam2bw_ATAC=${bam2bw_dir}/genomic_features/atac.bw
bam2bw_CTCF=${bam2bw_dir}/genomic_features/ctcf_log2fc.bw
zero_dir=${path}/data/mm10/${celltype}_zeroCTCF
zero_ATAC=${zero_dir}/genomic_features/atac.bw
zero_CTCF=${zero_dir}/genomic_features/ctcf_log2fc.bw
if [[ -f ${bam2bw_ATAC} ]] && [[ -f ${bam2bw_CTCF} ]]; then
echo ${celltype_name}
cat ${path}/hic_pred/marker_regions.txt | awk 'BEGIN{FS=OFS="\t"}$1!="chr"' | while read id
do
chr=`echo $id | awk '{FS="\t"}{print $1}'`
start=`echo $id | awk '{FS="\t"}{print $7}'`
echo ${chr}_${start}
python ~/softwares/C.Origami/src/corigami/inference/prediction.py --out "${path}/hic_pred/outputs/marker_regions" \
	--chr "${chr}" --celltype "${celltype}_bam2bw" --start ${start} --model "${path}/model_weights/corigami_base.ckpt" \
	--seq "${path}/data/mm10/dna_sequence" --ctcf ${bam2bw_CTCF} --atac ${bam2bw_ATAC}
sleep 1
done
else
echo "check ${celltype_name} CTCF files!"
sleep 3
fi
if [[ -f ${zero_ATAC} ]] && [[ -f ${zero_CTCF} ]]; then
echo ${celltype_name}
cat ${path}/hic_pred/marker_regions.txt | awk 'BEGIN{FS=OFS="\t"}$1!="chr"' | while read id
do
chr=`echo $id | awk '{FS="\t"}{print $1}'`
start=`echo $id | awk '{FS="\t"}{print $7}'`
echo ${chr}_${start}
python ~/softwares/C.Origami/src/corigami/inference/prediction.py --out "${path}/hic_pred/outputs/marker_regions" \
	--chr "${chr}" --celltype "${celltype}_zeroCTCF" --start ${start} --model "${path}/model_weights/corigami_base.ckpt" \
	--seq "${path}/data/mm10/dna_sequence" --ctcf ${zero_CTCF} --atac ${zero_ATAC}
sleep 1
done
else
echo "check ${celltype_name} zero files!"
sleep 3
fi
done

