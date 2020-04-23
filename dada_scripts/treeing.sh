#!/bin/bash

#$ -S /bin/bash
#$ -o treeing.stdout
#$ -e treeing.stderr 
#$ -l h_rt=48:00:00,h_vmem=6G
#$ -pe smp 10
#$ -binding linear:10
#$ -N treeing 

module load libgomp/4.4.7
export OMP_NUM_THREADS=10

clustalo="/data/project/metaamp/TOOLS/clustalo"
fasttree="/data/project/metaamp/TOOLS/FastTreeMP"

OTUfa=otus.fa
multAli=otus.multi.fa
outTree=tree.newick


date
echo "multiple alignment"
$clustalo -i $OTUfa -o $multAli --outfmt=fasta --threads=10 --force
date


echo "building tree"
$fasttree -nt -gamma -no2nd -fastest -spr 4 -log fasttree.log -quiet $multAli > $outTree
date

