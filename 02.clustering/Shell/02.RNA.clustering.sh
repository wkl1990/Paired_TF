#!/usr/bin/bash
#example: bash 02.RNA.clustering.sh -s XJ580 -a RNA -g mouse -e scRNA -p PairedTag/

tscc=/tscc/projects

while getopts ':s:ca:g:e:p:t:r:h' opt; do
  case "$opt" in
    s)
      sampleID="$OPTARG"
      echo "Processing option 's' with '${OPTARG}' argument"
      ;;
    c)
      combined="true"
      echo "Processing option 'c' with '${OPTARG}' argument"
      ;;
    a)
      assay="$OPTARG"
      echo "Processing option 'a' with '${OPTARG}' argument"
      ;;
    g)
      genome="$OPTARG"
      echo "Processing option 'g' with '${OPTARG}' argument"
      ;;
    e)
      environment="$OPTARG"
      echo "Processing option 'e' with '${OPTARG}' argument"
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
    h)
      echo "Usage: $(basename $0) [-s sampleID] [-c] [-g genome] [-m mode] [-p path]"
      exit 0
      ;;
    :)
      echo -e "option requires an argument.\nUsage: $(basename $0) [-s sampleID] [-c] [-a assay] [-g genome] [-e environment] [-p path]"
      exit 1
      ;;
    ?)
      echo -e "Invalid command option.\nUsage: $(basename $0) [-s sampleID] [-c] [-a assay] [-g genome] [-e environment] [-p path]"
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

if [ ! -d "21.seurat/cluster" ]; then
  mkdir -p 21.seurat/cluster/pca/
  mkdir -p 21.seurat/cluster/umap/
  mkdir -p 21.seurat/cluster/knn/
  mkdir -p 21.seurat/cluster/resolution/
  mkdir -p 21.seurat/cluster/silhouette/
  mkdir -p 21.seurat/cluster/cluster/
  mkdir -p 21.seurat/cluster/featurePlot/
  mkdir -p 21.seurat/cluster/marker/
  mkdir -p 21.seurat/cluster/log/
fi

# set environment
source ${tscc}/miniconda3/etc/profile.d/conda.sh
conda activate ${environment}

echo "STEP A: PCA"
if [[ $combined == "true" ]]; then
Rscript ${tscc}/scripts/R/02c.RNA.pc_estimate.R -i ${path}/21.seurat/preprocess/combined/${sampleID}.combined.rds -a ${assay} -c 100 -o ${path}/21.seurat/cluster/pca/${sampleID} &> ${path}/21.seurat/cluster/log/01.${sampleID}.pca.log
else
Rscript ${tscc}/scripts/R/02c.RNA.pc_estimate.R -i ${path}/21.seurat/preprocess/doublet/${sampleID}.doublet.rds -a ${assay} -c 100 -o ${path}/21.seurat/cluster/pca/${sampleID} &> ${path}/21.seurat/cluster/log/01.${sampleID}.pca.log
fi
Rscript ${tscc}/scripts/R/02c.RNA.findpc_auto.R -i ${path}/21.seurat/cluster/pca/${sampleID}.pca.rds -a ${assay} -c 100 -t -o ${path}/21.seurat/cluster/pca/${sampleID} &> ${path}/21.seurat/cluster/log/01.${sampleID}.autopc.log

echo "STEP B: UMAP"
npc=`cat ${path}/21.seurat/cluster/pca/${sampleID}.pc_all.csv | grep PCs | awk 'BEGIN{FS=","}{print $3}'`
Rscript ${tscc}/scripts/R/02d.RNA.umap.R -i ${path}/21.seurat/cluster/pca/${sampleID}.pca.rds -a ${assay} -c ${npc} -u raw -o ${path}/21.seurat/cluster/umap/${sampleID} &> ${path}/21.seurat/cluster/log/02.${sampleID}.umap.log

echo "STEP C: knn"
Rscript ${tscc}/scripts/R/02e.RNA.knn.R -i ${path}/21.seurat/cluster/umap/${sampleID}.umap.rds -r pca -a ${assay} -c 21 -o ${path}/21.seurat/cluster/knn/${sampleID} &> ${path}/21.seurat/cluster/log/03.${sampleID}.knn.log

echo "STEP D: consensus and silhouette plot"
if [[ -v parallel ]]; then
  bash ${tscc}/scripts/Shell/02f.consensusLeiden.sh ${path}/21.seurat/cluster/knn/${sampleID}.knn.mmtx ${path}/21.seurat/cluster/resolution/${sampleID} ${tscc} parallel &> ${path}/21.seurat/cluster/log/04.${sampleID}.resolution.log
else
  bash ${tscc}/scripts/Shell/02f.consensusLeiden.sh ${path}/21.seurat/cluster/knn/${sampleID}.knn.mmtx ${path}/21.seurat/cluster/resolution/${sampleID} ${tscc} &> ${path}/21.seurat/cluster/log/04.${sampleID}.resolution.log
fi
conda activate ${environment}

if [[ -v parallel ]]; then
  Rscript ${tscc}/scripts/R/02f.silhouette_resolution.R -i ${path}/21.seurat/cluster/knn/${sampleID}.knn.rds -m pca -a ${assay} -t -r -o ${path}/21.seurat/cluster/silhouette/${sampleID} &> ${path}/21.seurat/cluster/log/04.${sampleID}.silhouette.log
else
  Rscript ${tscc}/scripts/R/02f.silhouette_resolution.R -i ${path}/21.seurat/cluster/knn/${sampleID}.knn.rds -m pca -a ${assay} -t -o ${path}/21.seurat/cluster/silhouette/${sampleID} &> ${path}/21.seurat/cluster/log/04.${sampleID}.silhouette.log
fi

echo "STEP E: clustering"
if [[ -f ${path}/21.seurat/cluster/silhouette/${sampleID}.resolution_select.txt ]]; then
  resolution=`cat ${path}/21.seurat/cluster/silhouette/${sampleID}.resolution_select.txt`
elif [[ -f ${path}/21.seurat/cluster/resolution/${sampleID}.cluster.consensus.pdf ]]; then
  echo "consensus exist!"
  resolution=1.0
elif [[ -f ${path}/21.seurat/cluster/knn/${sampleID}.knn.rds ]]; then
  echo "consensus not exist!"
  resolution=1.0
else
  echo "knn not exist!"
fi

echo $resolution
Rscript ${tscc}/scripts/R/02g.RNA.clustering.R -i ${path}/21.seurat/cluster/knn/${sampleID}.knn.rds -r ${resolution} -a ${assay} -m 4 -t -o ${path}/21.seurat/cluster/cluster/${sampleID} &> ${path}/21.seurat/cluster/log/05.${sampleID}.clustering.log

echo "STEP F: plot feature"
Rscript ${tscc}/scripts/R/02i.RNA.brainfeature.R -i ${path}/21.seurat/cluster/cluster/${sampleID}.cluster.rds -a RNA -g ${genome} -s orig.ident -l Ctcf,Polr2a -o ${path}/21.seurat/cluster/featurePlot/${sampleID} &> ${path}/21.seurat/cluster/log/06.${sampleID}.featurePlot.log

echo "STEP G: pairwise marker"
Rscript ${tscc}/scripts/R/02j.RNA.pairmarker.R -i ${path}/21.seurat/cluster/cluster/${sampleID}.cluster.rds -a RNA -o ${path}/21.seurat/cluster/marker/${sampleID} &> ${path}/21.seurat/cluster/log/07.${sampleID}.pairwise.marker.log

echo "STEP H: all marker"
Rscript ${tscc}/scripts/R/02k.RNA.allmarker.R -i ${path}/21.seurat/cluster/cluster/${sampleID}.cluster.rds -a RNA -o ${path}/21.seurat/cluster/marker/${sampleID} &> ${path}/21.seurat/cluster/log/08.${sampleID}.all.marker.log

echo "Job is done!"

