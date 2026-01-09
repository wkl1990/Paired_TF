#!/usr/bin/bash
#example: bash 02.RNA.QC.sh -s XJ580 -g mouse -e scRNA -p PairedTag/

tscc=/tscc/projects
nfeature=500
mito=10

while getopts ':s:g:f:u:m:d:e:p:t:o:r:h' opt; do
  case "$opt" in
    s)
      sampleID="$OPTARG"
      echo "Processing option 's' with '${OPTARG}' argument"
      ;;
    g)
      genome="$OPTARG"
      echo "Processing option 'g' with '${OPTARG}' argument"
      ;;
    f)
      nfeature="$OPTARG"
      echo "Processing option 'f' with '${OPTARG}' argument"
      ;;
    u)
      ufeature="$OPTARG"
      echo "Processing option 'u' with '${OPTARG}' argument"
      ;;
    m)
      mito="$OPTARG"
      echo "Processing option 'm' with '${OPTARG}' argument"
      ;;
    d)
      dourate="$OPTARG"
      echo "Processing option 'd' with '${OPTARG}' argument"
      ;;
    e)
      environment="$OPTARG"
      echo "Processing option 'g' with '${OPTARG}' argument"
      ;;
    p)
      path="$OPTARG"
      echo "Processing option 'p' with '${OPTARG}' argument"
      ;;
    t)
      tscc="$OPTARG"
      echo "Processing option 't' with '${OPTARG}' argument"
      ;;
    r)
      parallel="$OPTARG"
      echo "Processing option 'r' with '${OPTARG}' argument"
      ;;
    o)
      out="$OPTARG"
      echo "Processing option 'o' with '${OPTARG}' argument"
      ;;
    h)
      echo "Usage: $(basename $0) [-s sampleID] [-g genome] [-f nfeature] [-m mode] [-p path]"
      exit 0
      ;;
    :)
      echo -e "option requires an argument.\nUsage: $(basename $0) [-s sampleID] [-g genome] [-f nfeature] [-e environment] [-p path]"
      exit 1
      ;;
    ?)
      echo -e "Invalid command option.\nUsage: $(basename $0) [-s sampleID] [-g genome] [-f nfeature] [-e environment] [-p path]"
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

raw_dir=${path}/21.seurat/preprocess/raw/
qc_dir=${path}/21.seurat/preprocess/qc/${out}
pca_dir=${path}/21.seurat/preprocess/pca/${out}
umap_dir=${path}/21.seurat/preprocess/umap/${out}
knn_dir=${path}/21.seurat/preprocess/knn/${out}
resolution_dir=${path}/21.seurat/preprocess/resolution/${out}
silhouette_dir=${path}/21.seurat/preprocess/silhouette/${out}
cluster_dir=${path}/21.seurat/preprocess/cluster/${out}
doublet_dir=${path}/21.seurat/preprocess/doublet/${out}
log_dir=${path}/21.seurat/preprocess/log/${out}

if [ ! -d "${path}/21.seurat/preprocess" ] || [[ $out != "" ]]; then
  mkdir -p ${raw_dir}/
  mkdir -p ${qc_dir}/
  mkdir -p ${pca_dir}/
  mkdir -p ${umap_dir}/
  mkdir -p ${knn_dir}/
  mkdir -p ${resolution_dir}/
  mkdir -p ${silhouette_dir}/
  mkdir -p ${cluster_dir}/
  mkdir -p ${doublet_dir}/
  mkdir -p ${log_dir}/
fi

if [[ "${genome}" == "mouse" ]]; then
  genome_version="mm10"
elif [[ "${genome}" == "human" ]]; then
  genome_version="hg38"
else
  echo -e "Error: ${genome} not found!"
  exit 1
fi

# set environment
source ${tscc}/miniconda3/etc/profile.d/conda.sh
conda deactivate
conda activate ${environment}

echo $tscc
which python
which Rscript
echo $PATH

echo "STEP A: Load data"
Rscript ${tscc}/scripts/R/02a.RNA.load_data.R -i ${path}/04.matrices/${sampleID}_${genome_version}_sorted_rmdup_mtx2 -n ${sampleID} -o ${raw_dir}/${sampleID} &> ${log_dir}/01.${sampleID}.load_data.log

echo "STEP B: QC"
if [[ -v ufeature ]]; then
  Rscript ${tscc}/scripts/R/02b.RNA.qc.R -i ${raw_dir}/${sampleID}.raw.rds -g ${genome} -n ${nfeature} -u ${ufeature} -m ${mito} -c 1 -t -o ${qc_dir}/${sampleID} &> ${log_dir}/02.${sampleID}.qc.log
else 
  Rscript ${tscc}/scripts/R/02b.RNA.qc.R -i ${raw_dir}/${sampleID}.raw.rds -g ${genome} -n ${nfeature} -m ${mito} -c 1 -t -o ${qc_dir}/${sampleID} &> ${log_dir}/02.${sampleID}.qc.log
fi

echo "STEP C: PCA"
Rscript ${tscc}/scripts/R/02c.RNA.pc_estimate.R -i ${qc_dir}/${sampleID}.qc.rds -a RNA -n 2000 -c 100 -o ${pca_dir}/${sampleID} &> ${log_dir}/03.${sampleID}.pca.log
Rscript ${tscc}/scripts/R/02c.RNA.findpc_auto.R -i ${pca_dir}/${sampleID}.pca.rds -a RNA -c 100 -t -o ${pca_dir}/${sampleID} &> ${log_dir}/03.${sampleID}.autopc.log

echo "STEP D: UMAP"
npc=`cat ${pca_dir}/${sampleID}.pc_all.csv | grep PCs | awk 'BEGIN{FS=","}{print $3}'`
echo $npc
Rscript ${tscc}/scripts/R/02d.RNA.umap.R -i ${pca_dir}/${sampleID}.pca.rds -a RNA -c $npc -u raw -o ${umap_dir}/${sampleID} &> ${log_dir}/04.${sampleID}.umap.log

echo "STEP E: knn"
npc=`cat ${pca_dir}/${sampleID}.pc_all.csv | grep PCs | awk 'BEGIN{FS=","}{print $3}'`
Rscript ${tscc}/scripts/02e.RNA.knn.R -i ${umap_dir}/${sampleID}.umap.rds -r pca -a RNA -c $npc -o ${knn_dir}/${sampleID} &> ${log_dir}/05.${sampleID}.knn.log

echo "STEP F: consensus and silhouette plot"
if [[ -v parallel ]]; then
  bash ${tscc}/scripts/Shell/02f.consensusLeiden.sh ${knn_dir}/${sampleID}.knn.mmtx ${resolution_dir}/${sampleID} ${tscc} parallel &> ${log_dir}/06.${sampleID}.resolution.log
else
  bash ${tscc}/scripts/Shell/02f.consensusLeiden.sh ${knn_dir}/${sampleID}.knn.mmtx ${resolution_dir}/${sampleID} ${tscc} &> ${log_dir}/06.${sampleID}.resolution.log
fi

if [[ -v parallel ]]; then
  Rscript ${tscc}/scripts/R/02f.silhouette_resolution.R -i ${knn_dir}/${sampleID}.knn.rds -m pca -a RNA -t -r -o ${silhouette_dir}/${sampleID} &> ${log_dir}/06.${sampleID}.silhouette.log
else
  Rscript ${tscc}/scripts/R/02f.silhouette_resolution.R -i ${knn_dir}/${sampleID}.knn.rds -m pca -a RNA -t -o ${silhouette_dir}/${sampleID} &> ${log_dir}/06.${sampleID}.silhouette.log
fi

echo "STEP G: clustering"
resolution=`cat ${silhouette_dir}/${sampleID}.resolution_select.txt`
echo $resolution
Rscript ${tscc}/scripts/R/02g.RNA.clustering.R -i ${knn_dir}/${sampleID}.knn.rds -r $resolution -a RNA -m 4 -o ${cluster_dir}/${sampleID} &> ${log_dir}/07.${sampleID}.clustering.log

echo "STEP H: DoubletFinder"
npc=`cat ${pca_dir}/${sampleID}.pc_all.csv | grep PCs | awk 'BEGIN{FS=","}{print $3}'`
if [[ -v dourate ]]; then
  Rscript ${tscc}/scripts/R/02h.RNA.doublet.R -i ${cluster_dir}/${sampleID}.cluster.rds -c $npc -d ${dourate} -m 1 -f -o ${doublet_dir}/${sampleID} &> ${log_dir}/08.${sampleID}.doublet.log
else 
  Rscript ${tscc}/scripts/R/02h.RNA.doublet.R -i ${cluster_dir}/${sampleID}.cluster.rds -c $npc -m 1 -f -o ${doublet_dir}/${sampleID} &> ${log_dir}/08.${sampleID}.doublet.log
fi

conda deactivate
echo "Job is done!"

