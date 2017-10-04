#!/bin/bash
#PBS -P xf3
#PBS -q express
#PBS -l walltime=24:00:00,mem=120GB,ncpus=16
#PBS -l jobfs=500GB

set -vx

#define some variables at the start
#give this job a name
name='Hu_FAH05731_SQK-LSK108_albacore202'
###
INPUT=/short/xf3/yh7166/MinION_data/20170617_FAH05731/raw_data
OUTPUT=/short/xf3/yh7166/MinION_data/20170617_FAH05731/basecalled_data
threads='16'
mem_size='120G'

#for basecalling change here
flowcellID="FAH05731"
kitID="SQK-LSK108"

#make the output folder
mkdir -p ${OUTPUT}

##this is a script to do basecalling and basic QC of your reads

#now move everything to the node so we can get started
cd $PBS_JOBFS
mkdir workspace
mkdir albacore_output
cp ${INPUT}/*tar.gz workspace/.

#go ahead with unziping and basecalling
cd workspace
for x in *.tar.gz
do
tar -xopf ${x}
done

mkdir fast5s
#check if this makes sense
mv */fast5 fast5s/.

#modules to load for basecalling
module load albacore/2.0.2

# basecall with albacore
read_fast5_basecaller.py -i ./fast5s -t $threads -s $PBS_JOBFS/albacore_output -f $flowcellID -k $kitID -r -n 0 -q 99999999999999999 -o fastq,fast5

#now pull out fastq and summary files before zipping up stuff

cd $PBS_JOBFS
mkdir ${name}
cp albacore_output/sequencing_summary.txt ${name}/.
cat albacore_output/*.fastq > ${name}/${name}.fastq


#remove original gz files and zip up stuff
rm -r workspace
tar -cvzf ${name}.tar.gz ${name}
rm -r ${name}

#now move everything from this step down already
mv ${name}.tar.gz ${OUTPUT}/.



#cd ${PBS_JOBFS}
#mv * ${OUTPUT}/.
