# What is dadasnake?
Dadasnake is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow to process amplicon sequencing data, from raw fastq-files to taxonomically assigned OTU tables, based on the [DADA2](http://benjjneb.github.io/dada2/) method. It is designed to run on a computing cluster using a single conda environment in multiple jobs triggered by Snakemake.

![overview](https://github.com/a-h-b/dadasnake/blob/master/documentation/pipeline.png)

## How to run dadasnake
To run the dadasnake, you need a config file and a sample table, plus data. The config file (in yaml format) is read by Snakemake to determine the inputs, steps, arguments and outputs. The sample table (tab-separated text) gives sample names and file names in the simplest case, with column headers named library and r1_file (and r2_file for paired-end data sets) - its path has to be mentioned in the config file. You can add columns labeled run and sample to indicate libraries that should be combined into one final column and different sequencing runs. All raw data (usually fastq files) need to be in the same directory (which has to be given in the config file). 
Once raw data, config file and sample file are present, the workflow can be started from the dadasnake directory by the snakemake command:
```
snakemake -s Snakefile --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda
```
If you're using a computing cluster, add your cluster's submission command and the number of jobs you want to maximally run at the same time, e.g.:
```
snakemake -j 50 -s Snakefile --cluster "qsub -l h_rt={params.runtime},h_vmem={params.mem} -pe smp {threads} -cwd" --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda
```
This will submit most steps as their own job to your cluster's queue.
If you want to share the conda installation with colleagues, use the `--conda-prefix` argument of Snakemake
```
snakemake -j 50 -s Snakefile --cluster "qsub -l h_rt={params.runtime},h_vmem={params.mem} -pe smp {threads} -cwd" --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda --conda-prefix /PATH/TO/YOUR/COMMON/CONDA/DIRECTORY
```
Depending on your dataset and settings, and your cluster's queue, the workflow will take a few minutes to days to finish. It can be helpful to run the Snakemake command in tmux. An example for this is provided in the dadasnake script in this repository.

## What does the dadasnake do?
* primer removal - using cutadapt
* quality filtering and trimming - using DADA2
* error estimation & denoising - using DADA2
* paired-ends assembly - using DADA2
* OTU table generation - using DADA2
* chimera removal - using DADA2
* taxonomic classification - using mothur and/or DECIPHER (& ITS detection - using ITSx & blastn)
* length check - in R
* treeing - using clustal omega and fasttree
* hand-off in biom-format, as R object, as R phyloseq object, and as fasta and tab-separated tables
* keeping tabs on number of reads in each step

## The configuration
The config file must be in .yaml format. 
