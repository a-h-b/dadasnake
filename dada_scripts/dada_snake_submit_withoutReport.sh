#! /bin/bash

#$ -cwd
#$ -N dadaSnake
#$ -l h_rt=172:00:00,h_vmem=8G
#$ -binding linear:1
#$ -e dadaSnake.stderr
#$ -o dadaSnake.stdout

CONFIGFILE=$1

module load miniconda/3/4.5.12-1
conda activate /data/project/metaamp/TOOLS/dada_pipe/dada_env >/dev/null 2>&1

snakemake -j 50 -s /data/project/metaamp/TOOLS/dada_pipe/Snakefile --cluster "qsub -l h_rt={params.runtime},h_vmem={params.mem} -pe smp {threads} -cwd" --configfile $CONFIGFILE

conda deactivate >/dev/null 2>&1
