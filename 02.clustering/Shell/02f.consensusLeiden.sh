#!/usr/bash

#tscc=/projects

#if [ ! -d "resolution" ]; then
#	mkdir resolution
#fi

input=${1} 
output=${2}
tscc=${3}
#parallel=${4}


echo "Step1: consensus analysis!"

if [ -n "$4" ]; then
	for r in `seq 0.1 0.1 2`;
	do echo "Resolution: "${r}
	#python $path2script/snapATAC.leiden.py -i hba.cortex.knn.mmtx -r ${r} -o refineCluster/hba.cortex.${r}
	#Rscript $path2script/snapATAC.refineCluster.R -i hba.cortex.cluster.RData -p refineCluster/hba.cortex.${r}.partition.txt -o refineCluster/hba.cortex.${r}
	python ${tscc}/scripts/Python/02f.consensusAnalysis.py -i ${input} -r ${r} -n 100 -o ${output}.${r} &
	sleep 10
	done
	wait
else
	for r in `seq 0.1 0.1 2`;
	do echo "Resolution: "${r}
	#python $path2script/snapATAC.leiden.py -i hba.cortex.knn.mmtx -r ${r} -o refineCluster/hba.cortex.${r}
	#Rscript $path2script/snapATAC.refineCluster.R -i hba.cortex.cluster.RData -p refineCluster/hba.cortex.${r}.partition.txt -o refineCluster/hba.cortex.${r}
	python ${tscc}/scripts/Python/02f.consensusAnalysis.py -i ${input} -r ${r} -n 100 -o ${output}.${r} 
	done
fi


echo "Step2: generate consensus statistics!"
if [ ! -f ${output}.consensus.stat.txt ]; then
	cat ${output}.*.stat.txt > ${output}.consensus.stat.txt
	echo "Step3: consensus plot!"
	Rscript ${tscc}/scripts/R/02f.consensusPlot.R -i ${output} -t ${output}
else
	echo "Error: consensus statistics exist!"
fi
echo "Job is done!"



