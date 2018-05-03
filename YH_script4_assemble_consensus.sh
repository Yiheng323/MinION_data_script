#!/bin/bash
##SBATCH--job-name=consensus
#SBATCH--output=cons
#SBATCH--nodes=1
#SBATCH--ntasks-per-node=8
#SBATCH--mem=32G#mb


###Changethefollowingvariablesaccordingtoyoursettings.
RACON=/home/sonnguyen/sw/racon/# Your installed racon path
INPUT=/home/merged.fastq # Merged fastq file,to merge fastq files into a single fastq file, go to the directory of the fastq files, and type "cat *.fastq > merge.fastq"
OUTPUT=/home/output#Output_folder
###
minimap2 -t 8 -x ava-ont ${INPUT} ${INPUT} > ${OUTPUT}/reads.paf

miniasm -f ${INPUT} ${OUTPUT}/reads.paf > ${OUTPUT}/raw_contigs.gfa

awk'$1~/S/{print">"$2"\n"$3}' ${OUTPUT}/raw_contigs.gfa > ${OUTPUT}/raw_contigs.gfa

minimap2 ${OUTPUT}/raw_contigs.gfa ${INPUT} > ${OUTPUT}/mapping.paf

${RACON}/bin/racon ${INPUT} ${OUTPUT}/mapping.paf ${OUTPUT}/raw_contigs.gfa ${OUTPUT}/consensus_0.fasta#

##raconloop
for j in 0 1 2 3 4 5 6 7 8 9 10
  do
    minimap2 ${OUTPUT}/consensus_${j}.fasta ${INPUT} | ${RACON}/bin/racon ${INPUT} - ${OUTPUT}/consensus_${j}.fasta ${OUTPUT}/consensus_$((j+1)).fasta
  done

done
