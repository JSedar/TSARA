#!/bin/bash
#mkdir -p local
#conda activate nest
export LC_ALL="en_GB.UTF-8"
workDir=/home/samara/NeST/Tsara
outDir=/home/samara/Bureau/Resultats
myFolder=Tsara_$(date "+%Y_%m_%d_%s")
time python3 /home/samara/NeST/nest.py \
	-i ${workDir}/FichierFastq_3Aout \
	-a ${workDir}/adapters.fa \
	-r ${workDir}/ref_fa.fasta \
	-b ${workDir}/ref_bed.bed \
	-o ${outDir}/${myFolder} \
	--varofint ${workDir}/knownSNP1.csv

time python3 /home/samara/NeST/Tsara/Depth_HRP.py ${outDir}/${myFolder}