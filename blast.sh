
#!/bin/bash
#$ -M yiheng.hu@anu.edu.au
#$ -m a
#$ -cwd
#$ -V
#$ -j y
#$ -pe threads 16
#$ -l virtual_free=3.0g,h_vmem=3.0g
#$ -N _blast
set -vx
#uses a local nt database downloaded 20170802 to blast all contigs against this database
WORKDIR=~/data/20170617_FAH05731/ncbiblast
for fasta in ~/data/20170617_FAH05731/ncbiblast/*.fasta; do blastn -query $fasta -db nt -evalue 0.01 -outfmt '6 qseqid sseqid evalue bitscore length pident nident sgi sacc staxids sscinames scomnames sskingdoms' -show_gis -num_threads 16 | sort -k1,1 -k4,4nr | sort -u -k1,1 --merge > $fasta.20170807_NCBI.nt.ncbiblast_output; done
