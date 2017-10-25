#######################################################
#YH_script_1:from albacore output to stats information#
#######################################################

### the default operation cluster is GDU

### Before start:
### 1. install anaconda:
###    download (wget) Anaconda_XXX.sh
###    bash Anaconda_XXX.sh
###    vim .bashrc, insert "unset PYTHONPATH" to unset the default python path of GDU, to make anaconda avalaible
### 2. set up bioconda in GDU:
###    conda config --add channels defaults
###    conda config --add channels conda-forge
###    conda config --add channels bioconda
### 3. install NanoLyse, porechop, bbmap:
###    conda install NanoLyse
###    conda install bbmap
###    git clone https://github.com/rrwick/Porechop.git
###    cd Porechop
###    python3 setup.py install
###    porechop -h
### 4. make rg database containing the species of interests, in this study the rg are made from wheat, Pst, Zymo strain W332, Pyre and stago reference genomes (all used unmasked genomes)
###    (1). cd ~/bio/
###    (2). download wheat genome from EnsemblPlants: ftp://ftp.ensemblgenomes.org/pub/plants/release-37/fasta/triticum_aestivum/dna/Triticum_aestivum.TGACv1.dna.nonchromosomal.fa.gz
###    (2). download Pst, Zymo, Pyre, Stago reference sequence from MycoCosm:
###         Pst: https://genome.jgi.doe.gov/Pucstr1/download/Pucstr1_GeneCatalog_proteins_20170922_promoters_1k.fa.gz
###         Zymo: https://genome.jgi.doe.gov/Mycgr3/download/Mycosphaerella_graminicola.fasta.gz
###         Pyre: https://genome.jgi.doe.gov/Pyrtr1/download/Pyrenophora_tritici_repentis_unmasked_assembly.fasta.gz
###         Stago: https://genome.jgi.doe.gov/Stano2/download/Stagonospora_nodorum.fasta
###    (3). gunzip all of the .gz files
###    (4). change the header of each reference genome to be their names:
###         awk '/^>/{print ">Wheat_" ++i; next}{print}' < Triticum_aestivum.TGACv1.dna.nonchromosomal.fa > Wheat.fasta
###         awk '/^>/{print ">Stago_" ++i; next}{print}' < stagonospora_nodorum.fasta > Stago.fasta
###         awk '/^>/{print ">Zymo_" ++i; next}{print}' < Mycosphaerella_graminicola.fasta > Zymo.fasta
###         awk '/^>/{print ">Pyre_" ++i; next}{print}' < Pyrenophora_tritici_repentis_unmasked_assembly.fasta > Pyre.fasta
###         awk '/^>/{print ">Pst_" ++i; next}{print}' < Pucstr1_GeneCatalog_proteins_20170922_promoters_1k.fa > Pst.fasta
###    (5). cat Wheat.fasta Stago.fasta Zymo.fasta Pyre.fasta Pst.fasta > rg.fasta
###    (6). makeblastdb -in rg.fasta -dbtype nucl -title rg -out rg
###           Building a new DB, current time: 10/24/2017 17:03:11
###           New DB name:   rg
###           New DB title:  rg
###           Sequence type: Nucleotide
###           Keep Linkouts: T
###           Keep MBits: T
###           Maximum file size: 1000000000B
###           Adding sequences from FASTA; added 756517 sequences in 236.097 seconds.




#!/bin/bash
#$ -M yiheng.hu@anu.edu.au
#$ -m a
#$ -cwd
#$ -V
#$ -j y
#$ -pe threads 12
#$ -l h_vmem=3g,virtual_free=2.9g
#$ -N pipe1
set -vx

BASEFOLDER='/home/yiheng/test' # make sure the input file is only a copy
NAME='Hu_FAH05731_albacore202'
DATE='20171025'
BARCODE='barcode10'

cd ${BASEFOLDER}/basecalled_data
tar -xvf ${NAME}.tar.gz
cd ${NAME}/workspace

### filter CDS

cat pass/${BARCODE}/*.fastq fail/${BARCODE}/*.fastq > ${NAME}.${BARCODE}.unlysed.fastq
gzip ${NAME}.${BARCODE}.unlysed.fastq
gunzip -c ${NAME}.${BARCODE}.unlysed.fastq.gz | NanoLyse | gzip > ${NAME}.${BARCODE}.fastq.gz
gunzip ${NAME}.${BARCODE}.fastq.gz
rm ${NAME}.${BARCODE}.unlysed.fastq.gz

### move to the workspace folder for further manipulation
cd ${BASEFOLDER}
mkdir -p workspace/${BARCODE}
mv ${BASEFOLDER}/basecalled_data/${NAME}/workspace/pass/${NAME}.${BARCODE}.fastq ${BASEFOLDER}/workspace/${BARCODE}
cd ${BASEFOLDER}/workspace/${BARCODE}

# do porechop to chop out adapter sequence

porechop -i ${BASEFOLDER}/workspace/${BARCODE}/${NAME}.${BARCODE}.fastq -o ${BASEFOLDER}/workspace/${BARCODE}/${NAME}.chopped.${BARCODE}.fastq --format fastq --middle_threshold 95

# convert the fastq to fasta

sed '/^@/!d;s//>/;N' ${BASEFOLDER}/workspace/${BARCODE}/${NAME}.chopped.${BARCODE}.fastq > ${BASEFOLDER}/workspace/${BARCODE}/${NAME}.chopped.${BARCODE}.fasta

# do blastn for the fasta file

blastn -query ${NAME}.chopped.${BARCODE}.fasta -db rg -evalue 0.01 -outfmt '6 qseqid sseqid evalue bitscore length pident nident sgi sacc staxids sscinames scomnames sskingdoms' -show_gis -num_threads 12 | sort -k1,1 -k4,4nr | sort -u -k1,1 --merge > ${NAME}_${BARCODE}_chopped.fasta.${DATE}.rgblast_output

### supplimentray information for format6:

### qseqid means Query Seq-id
### sseqid means Subject Seq-id
### evalue means Expect value
### bitscore means Bit score
### length means Alignment length
### pident means Percentage of identical matches
### nident means Number of identical matches
### sgi means Subject GI
### sacc means Subject accession
### staxids means Subject Taxonomy ID(s), separated by a ';'
### sscinames means Subject Scientific Name(s), separated by a ';'
### scomnames means Subject Common Name(s), separated by a ';'
### sskingdoms means Subject Super Kingdom(s), separated by a ';'

# cut the seqid from rgblast_output and write it into a list for separate the fasta hit files

cut -f 1 ${NAME}_${BARCODE}_chopped.fasta.${DATE}.rgblast_output > ${NAME}.rgblast.qseqid.${BARCODE}.txt

# separate the blast hit and nohit reads from the fasta file

filterbyname.sh in=${NAME}.chopped.${BARCODE}.fasta out=${NAME}.chopped.rghityes.${BARCODE}.fasta names=${NAME}.rgblast.qseqid.${BARCODE}.txt include=t
filterbyname.sh in=${NAME}.chopped.${BARCODE}.fasta out=${NAME}.chopped.rghitno.${BARCODE}.fasta names=${NAME}.rgblast.qseqid.${BARCODE}.txt include=f

# do the second blast for the rghitno fasta file

blastn -query ${NAME}.chopped.rghitno.${BARCODE}.fasta -db nt -evalue 0.01 -outfmt '6 qseqid sseqid evalue bitscore length pident nident sgi sacc staxids sscinames scomnames sskingdoms' -show_gis -num_threads 12 | sort -k1,1 -k4,4nr | sort -u -k1,1 --merge > ${NAME}_${BARCODE}_chopped.fasta.${DATE}.ntblast_output

# cut the seqid from the ntblast_output and write it into a list for separate the fasta hit files

cut -f 1 ${NAME}_${BARCODE}_chopped.fasta.${DATE}.ntblast_output > ${NAME}.ntblast.qseqid.${BARCODE}.txt

# separate the blast hit and nohit reads from the fasta file

filterbyname.sh in=${NAME}.chopped.rghitno.${BARCODE}.fasta out=${NAME}.chopped.rghitno.nthityes.${BARCODE}.fasta names=${NAME}.ntblast.qseqid.${BARCODE}.txt include=t
filterbyname.sh in=${NAME}.chopped.rghitno.${BARCODE}.fasta out=${NAME}.chopped.rghitno.nthitno.${BARCODE}.fasta names=${NAME}.ntblast.qseqid.${BARCODE}.txt include=f
