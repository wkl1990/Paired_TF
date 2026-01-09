#!/bin/bash
path=/tscc/projects


PATHTO=/tscc/softwares/ROSE/
PYTHONPATH=$PATHTO/lib
export PYTHONPATH
export PATH=$PATH:$PATHTO/bin

cat ${path}/21.seurat/cluster/RNAPII_meta.table.sampleID.csv | awk 'BEGIN{FS=","}NR>1{print $8}' | sort | uniq | while read name 
do
echo $name
celltype=${name#*_}
echo $celltype
nohup ROSE_main.py -g MM10 -i ${path}/RNAPII/peak_gff/${celltype}.gff -r ${path}/RNAPII/bam/${celltype}.bam -o ${path}/RNAPII/rose/${celltype} -s 12500 -t 0 2>&1 > ${path}/RNAPII/log/${celltype}.ROSE.s12500.t0.log &
sleep 3
loadavg
done
