![badge](https://zenodo.org/badge/208493040.svg)

![logo](https://github.com/a-h-b/dadasnake/blob/master/documentation/snake_all-trans_b.png)

dadasnake is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow to process amplicon sequencing data, from raw fastq-files to taxonomically assigned "OTU" tables, based on the [DADA2](http://benjjneb.github.io/dada2/) method. Running dadasnake could not be easier: it is called by a single command from the command line. With a human-readable configuration file and a simple sample table, its steps are adjustable to a wide array of input data and requirements. It is designed to run on a computing cluster using a single conda environment in multiple jobs triggered by Snakemake. dadasnake reports on intermediary steps and statistics in intuitive figures and tables. Final data output formats include biom format, phyloseq objects, and flexible text files or R data sets for easy integration in microbial ecology analysis scripts.

## Installing dadasnake
For dadasnake to work, you need [conda](https://www.anaconda.com/). 

1) Clone this repository to your disk:
```
git clone https://github.com/a-h-b/dadasnake.git
```
Change into the dadasnake directory:
```
cd dadasnake
```
At this point, you have all the scripts you need to run the workflow using snakemake, and you'd just need to get some data and databases (see point 8). If you want to use the **comfortable dadasnake wrapper**, follow the points 2-6. 

2) Decide how you want to run dadasnake, if you let it submit jobs to the cluster:
Only do one of the two:
* if you want to submit the process running snakemake to the cluster:
```
cp auxiliary_files/dadasnake_allSubmit dadasnake
chmod 755 dadasnake
```
* if you want to keep the process running snakemake on the frontend using tmux:
```
cp auxiliary_files/dadasnake_tmux dadasnake
chmod 755 dadasnake
```

If you don't submit jobs to the cluster, but want to run the whole workflow interactively, e.g. on a laptop, it doesn't matter which wrapper you use. Just copy one of them, as described above.

3) Adjust the file VARIABLE_CONFIG to your requirements (have a tab between the variable name and your setting):
* SNAKEMAKE_VIA_CONDA - set this to true, if you don't have snakemake in your path and want to install it via conda. Leave empty, if you don't need an additional snakemake.
* SNAKEMAKE_EXTRA_ARGUMENTS - if you want to pass additional arguments to snakemake, put them here (e.g. --latency-wait=320 for slower file systems). Leave empty usually. 
* LOADING_MODULES - insert a bash command to load modules, if you need them to run conda. Leave empty, if you don't need to load a module.
* SUBMIT_COMMAND - insert the bash command you'll usually use to submit a job to your cluster to run on a single cpu for a few days. You only need this, if you want to have the snakemake top instance running in a submitted job. You alternatively have the option to run it on the frontend via tmux. Leave empty, if you want to use this frontend version and have [tmux](https://github.com/tmux/tmux/wiki) installed. You don't need to set this, if you are wanting to run the workflow interactively / on a laptop.
* BIND_JOBS_TO_MAIN - if you use the option to run the snakemake top instance in a submitted job and need to bind the other jobs to the same node, you can set this option to true. See FAQ below for more details. You don't need to set this, if you are wanting to run the workflow interactively / on a laptop.
* NODENAME_VAR - if you use the BIND_JOBS_TO_MAIN option, you need to let dadasnake know, how to access the node name (e.g.SLURMD_NODENAME on slurm). You don't need to set this, if you are wanting to run the workflow interactively / on a laptop.
* SCHEDULER - insert the name of the scheduler you want to use (currently `slurm` or `uge`). This determines the cluster config given to snakemake, e.g. the cluster config file for slurm is config/slurm.config.yaml . Also check that the settings in this file is correct. If you have a different system, contact us ( https://github.com/a-h-b/dadasnake/issues ). You don't need to set this, if you are wanting to run the workflow interactively / on a laptop.
* MAX_THREADS - set this to the maximum number of cores you want to be using in a run. If you don't set this, the default will be 50. Users can override this setting at runtime.
* NORMAL_MEM_EACH - set the size of the RAM of one core of your normal copute nodes (e.g. 8G). If you're not planning to use dadasnake to submit to a cluster, you don't need to set this. 
* BIGMEM_MEM_EACH - set the size of the RAM of one core of your bigmem (or highmem) compute nodes. If you're not planning to use dadasnake to submit to a cluster or don't have separate bigmem nodes, you don't need to set this.
* BIGMEM_CORES - set this to the maximum number of bigmem cores you want to require for a task. Set to 0, if you don't have separate bigmem nodes. You don't need to set this, if you're not planning to use dadasnake to submit to a cluster.
* LOCK_SETTINGS - set this to true, if you don't want users to choose numbers and sizes of compute nodes at run time. If you're not planning to use dadasnake to submit to a cluster, you don't need to set this. Setting LOCK_SETTINGS  makes the workflow slightly less flexible, as all large data sets will be run with the maximum number of bigmem nodes you set up here (see big_data settings below). On the other hand, it can be helpful, if you're setting up dadasnake for inexperienced users or have only one possible setting anyhow. If you're not locking, it's advised to set useful settings in the config/config.default.yaml file for normalMem, bigMem, and bigCores.


4) **optional, but highly recommended**: Install snakemake via conda:
If you want to use snakemake via conda (and you've set SNAKEMAKE_VIA_CONDA to true), install the environment, as [recommended by Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):
```
conda install -c conda-forge mamba
mkdir -p conda
mamba create --prefix $PWD/conda/snakemake_env
conda activate $PWD/conda/snakemake_env
mamba install -c conda-forge -c bioconda snakemake=6.9.1 mamba
conda deactivate
```
Alternatively, if the above does not work, you can install a fixed snakemake version without mamba like so:
```
conda env create -f workflow/envs/snakemake_env.yml --prefix $PWD/conda/snakemake_env
```
Dadasnake will run with Snakemake version >= 5.9.1 and hasn't been tested with any previous versions.

5) Set permissions / PATH:
Dadasnake is meant to be used by multiple users. Set the permissions accordingly. I'd suggest:
* to have read access for all files for the users plus 
* execution rights for the dadasnake file and the .sh scripts in the subfolder submit_scripts
* read, write and execution rights for the conda subfolder 
* Add the dadasnake directory to your path. 
* It can also be useful to make the VARIABLE_CONFIG file not-writable, because you will always need it. The same goes for config.default.yaml once you've set the paths to the databases you want to use (see below).

6) Initialize conda environments:
This run sets up the conda environments that will be usable by all users:
```
./dadasnake -i config/config.init.yaml 
```
This step will take several minutes. It will also create a folder with the name "dadasnake_initialized". You can safely remove it or keep it.
I strongly suggest to **remove one line from the activation script** after the installation, namely the one reading: `R CMD javareconf > /dev/null 2>&1 || true`, because you don't need this line later and if two users run this at the same time it can cause trouble. You can do this by running:
```
sed -i "s/R CMD javareconf/#R CMD javareconf/" conda/*/etc/conda/activate.d/activate-r-base.sh
```

7) **Optional** test run:
The test run does not need any databases. You should be able to start it by running 
```
./dadasnake -l -n "TESTRUN" -r config/config.test.yaml
```
If all goes well, dadasnake will run in the current session, load the conda environment, and make and fill a directory called testoutput. A completed run contains a file "workflow.done". 
If you don't want to see dadasnake's guts at this point, you can also run this with the -c or -f settings to submit to your cluster or start a tmux session (see How to run dadasnake below). 

8) Databases:
The dadasnake does not supply databases. I'd suggest to use the [SILVA database](https://www.arb-silva.de/no_cache/download/archive/current/Exports/) for 16S data and [UNITE](https://doi.org//10.15156/BIO/786336) for ITS. 
* dadasnake can use [mothur](https://www.mothur.org/) to do the classification, as it's faster and likely more accurate than the legacy DADA2 option. You need to format the database like for mothur ([see here](https://www.mothur.org/wiki/Taxonomy_outline)). 
* dadasnake can alternatively use the DADA2 implementation of the same classifier. You can find some databases maintained by Michael R. McLaren [here](https://zenodo.org/record/3986799). More information on the format is in the [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html).
* In addition to the bayesian classifier, dadasnake implements [DECIPHER](http://www2.decipher.codes/Documentation.html). You can find decipher [databases](http://www2.decipher.codes/Downloads.html) on the decipher website or build them yourself. 
* dadasnake can use [fungal traits](https://link.springer.com/article/10.1007/s13225-020-00466-2#data-availability) to assign traits to fungal genere. Download the latest table from [here](https://docs.google.com/spreadsheets/d/1cxImJWMYVTr6uIQXcTLwK1YNNzQvKJJifzzNpKCM6O0/edit#gid=492619054) - dadasnake has been tested with v1.2.
* You can also use dadasnake to blast and summarize results using [basta](https://github.com/timkahlke/BASTA). Have a look at the [NCBI's ftp](https://ftp.ncbi.nlm.nih.gov/blast/db/).
* To annotate fungal taxonomy with guilds via funguild, if you have suitable databases. 
* If you still have a tax4fun2 installation, you can also use it within dadasnake. The package and database were taken off github, so it's not part of defaut dadasnake anymore. 
* You can now use [picrust2](https://github.com/picrust/picrust2/wiki/) within dadasnake.
**You need to set the path to the databases of your choice in the config file.** By default, dadasnake looks for databases in the directory above where it was called. It makes sense to change this for your system in the config.default.yaml file upon installation, if all users access databases in the same place.

9) Fasttree:
dadasnake comes with fasttree for treeing, but if you have a decent number of sequences, it is likely to be relatively slow. If you have fasttreeMP, you can give the path to it in the config file.

## How to cite dadasnake
[Christina Weißbecker, Beatrix Schnabel, Anna Heintz-Buschart, Dadasnake, a Snakemake implementation of DADA2 to process amplicon sequencing data for microbial ecology, GigaScience, Volume 9, Issue 12, December 2020, giaa135](https://doi.org/10.1093/gigascience/giaa135). Please also cite [DADA2](https://www.nature.com/articles/nmeth.3869): Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016), and any other tools you use within dadasnake, e.g. [mothur](https://mothur.org/wiki/frequently_asked_questions/#how-do-i-cite-mothur-how-about-the-individual-functions), [DECIPHER](http://www2.decipher.codes/Citation.html), [ITSx](https://microbiology.se/software/itsx), [Fasttree](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490), [VSEARCH](https://doi.org/10.7717/peerj.2584), [FUNGuild](https://www.sciencedirect.com/science/article/abs/pii/S1754504815000847), [PICRUSt2](https://doi.org/10.1038/s41587-020-0548-6), [BASTA](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13095), [tax4fun2](https://environmentalmicrobiome.biomedcentral.com/articles/10.1186/s40793-020-00358-7).

##  

![overview](https://github.com/a-h-b/dadasnake/blob/master/documentation/pipeline.png)

## How to run dadasnake
To run the dadasnake, you need a config file and a sample table, plus data: 
* The config file (in yaml format) is read by Snakemake to determine the inputs, steps, arguments and outputs. 
* The sample table (tab-separated text) always gives sample names and file names, with column headers named library and r1_file (and r2_file for paired-end data sets). The path to the sample table has to be mentioned in the config file. You can add columns labeled `run` and `sample` to indicate libraries that should be combined into one final column and different sequencing runs (see the section about the sample table below). 
* All raw data (usually fastq files) need to be in one directory (which has to be given in the config file). 
* It is possible (and the best way to do this) to have one config file per run, which defines all settings that differ from the default config file.

### Using the dadasnake wrapper
As shown in the installation description above, dadasnake can be run in a single step, by calling dadasnake. Since most of the configuration is done via the config file, the options are very limited. You can either:
* -c run (submit to a cluster) dadasnake and make a report (-r), or
* -l run (in the current terminal) dadasnake and make a report (-r), or
* -f run (in a tmux session on the frontend) dadasnake *only available in the tmux installation* and make a report (-r), or
* just make a report (-r), or 
* run a dryrun (-d), or 
* unlock a working directory, if a run was killed (-u)
* initialize the conda environmnets only (-i) - you should only need this during the installation. 
It is strongly recommended to **first run a dryrun on a new configuration**, which will tell you within a few seconds and without submission to a cluster whether your chosen steps work together, the input files are where you want them, and your sample file is formatted correctly. In all cases you need the config file as the last argument. 
```
dadasnake -d -r config.yaml
```
You can also set the number of cpus to maximally run at the same time with -t. The defaults (1 for local/frontend runs and 50 for clusters) are reasonable for many settings and if you don't know what this means, you probably don't have to worry. But you may want to increase the numbers for larger datasets or bigger infrastructure, or decrease the numbers to match your environment's constraints.
You can add a name for your main job (-n NAME), e.g.:
```
dadasnake -c -n RUNNAME -r config.yaml
```
Note that spaces in RUNNAME are not allowed and dots will be replaced by underscores.

If you use the tmux version, you can see the tmux process running by typing `tmux ls`. You can also see the progress by checking the stdandard error file `tail RUNNAME_XXXXXXXXXX.stderr`.

Depending on your dataset and settings and your cluster's scheduler, the workflow will take a few minutes to days to finish. 

### Running snakemake manually
Once raw data, config file and sample file are present, the workflow can be started from the dadasnake directory by the snakemake command:
```
snakemake -s Snakefile --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda
```
If you're using a computing cluster, add your cluster's submission command and the number of jobs you want to maximally run at the same time, e.g.:
```
snakemake -j 50 -s Snakefile --cluster "qsub -l h_rt={resources.runtime},h_vmem=8G -pe smp {threads} -cwd" --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda 
```
This will submit most steps as their own job to your cluster's queue. The same can be achieved with a [cluster configuration](https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html#cluster-execution):
```
snakemake -j 50 -s Snakefile --cluster-config PATH/TO/SCHEDULER.config.yaml --cluster "{cluster.call} {cluster.runtime}{resources.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.threads}{threads} {cluster.partition}" --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda
```
If you want to share the conda installation with colleagues, use the `--conda-prefix` argument of Snakemake
```
snakemake -j 50 -s Snakefile --cluster-config PATH/TO/SCHEDULER.config.yaml --cluster "{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.threads}{threads} {cluster.partition}" --use-conda --conda-prefix /PATH/TO/YOUR/COMMON/CONDA/DIRECTORY
```
Depending on your dataset and settings, and your cluster's queue, the workflow will take a few minutes to days to finish.

## What does the dadasnake do?
* primer removal and removal of poly-G-tails - using cutadapt
* quality filtering and trimming - using DADA2
* optional downsampling of reads per sample - using [seqtk](https://github.com/lh3/seqtk) 
* error estimation & denoising - using DADA2, including Novaseq-enabled models
* paired-ends assembly - using DADA2
* "OTU" table generation (it contains ASVs, of course) - using DADA2
* chimera removal - using DADA2
* clustering of ASVs at a user-set similarity (these are called OTU now)
* taxonomic classification - using mothur and/or DECIPHER (& ITS detection - using ITSx & blastn + BASTA)
* functional annotation - using funguild, fungalTraits, picrust2 (or tax4fun2)
* length check - in R
* treeing - using clustal omega and fasttree
* hand-off in biom-format, as R object, as R phyloseq object, and as fasta and tab-separated tables
* keeping tabs on number of reads in each step, and read quality control - using fastqc & multiQC
You can control the settings for each step in a config file.

![steps](https://github.com/a-h-b/dadasnake/blob/master/documentation/steps.png)

## The samples table
Every samples table needs sample names (under header library) and file names (just the names, the path should be in the config file under header r1_file and potentially r2_file). Since DADA2 estimates run-specific errors, it can be helpful to give run IDs (under header run). If you have many (>500 samples), it is also useful to split them into runs for the analysis, as some of the most memory-intensive steps are done by run.  
If several fastq files should end up in the same column of the ASV/OTU table, you can indicate this by giving these libraries the same sample name (under header sample). Libraries from different runs are combined in the final ASV/OTU table (example 1). Libraries from the same run are combined after primer-processing (example 2).
Example 1:
![overview](https://github.com/a-h-b/dadasnake/blob/master/documentation/samples_ex1.png)
Example 2:
![overview](https://github.com/a-h-b/dadasnake/blob/master/documentation/samples_ex2.png)


## The configuration
The config file must be in .yaml format. The order within the yaml file does not matter, but the hierarchy has to be kept. Explanations can be found in the config-file in config/config.default.yaml .

**top-level parameters** | **sub-parameters** | **subsub-parameters** | **default value** | **possible values** | **used in stage** | **explanation** | **comments / recommendations**
---|---|---|---|---|---|---|---
raw_directory |  |  | "testdata"| any one path where you might have your raw data | all | directory with all raw data | you will usually have this somewhere in a project folder
sample_table |  |  | "testdata/samples.small.tsv" | any one location of your samples table | all | path to the samples table | the dadasnake will copy it to your output directory
outputdir |  |  | "dadasnake_output" | any path that you have permissions for writing to | all | directory where all the output will go | change this; a scratch-type place works best; each output directory can hold the results of one completed pipeline run only
paired||| true|true or false|primers and dada|do you want to use paired-end sequencing data?|if true, you have to give r1_file and r2_file in the samples table, if false only r1_file is read (if you want to use only R2 files from a paired-end sequencing run, put their name in the r1_file column)
tmp_dir |  |  | "tmp" | any path that you have permissions for writing to | all | directory for temporary, intermediate files that shouldn't be kept | keep this in a temporary place so you don't need to worry about removing its contents
big_data |  |  | false | a boolean | dada, taxonomy, post | whether to use big data settings | set to true, if you have extra high memory nodes and more than 1000 samples 
email |  |  | "" | "" or a valid email address | all | email address for mail notification | keep empty if you don't want emails. Check spelling, it's not tested.
do_primers |  |  | true | true or false | all | should primers be cut? | 
do_dada |  |  | true | true or false | all | should DADA2 be run? | 
do_taxonomy |  |  | true | true or false | all | should taxonomic classification be done? |
do_postprocessing |  |  | true | true or false | all | should some more steps be done (e.g. functional annotation) | 
hand_off|||||dada, taxonomy, postprocessing||settings deciding if additional formats should be given
&nbsp;|  biom||false|true or false|dada, taxonomy|whether a biome format output should be written|biome contains ASV table or ASV table and taxonomy (if taxonomy was run); biome table is never filtered
&nbsp;|  phyloseq||true|true or false|taxonomy, postprocessing|whether a phyloseq object (or two - for ASVs and OTUs) should be returned|contains ASV or OTU table and taxonomy and tree (if each was run; if tree is run on pruned OTU or ASV table, phyloseq object contains filtered dataset)
primers |  |  |  |  | primers |  | information on primers
  &nbsp;| fwd |  |  |  | primers |  | information on forward primer
  &nbsp;|  | sequence | GTGYCAGCMGCCGCGGTAA | any sequence of IUPAC DNA code | primers | sequence of forward primer |
 &nbsp;|  | name | 515F | anything | primers | name of forward primer | for your reference only
&nbsp;| rvs |  |  |  | primers |  | information on reverse primer
&nbsp;||    sequence| GGACTACNVGGGTWTCTAAT|any sequence of IUPAC DNA code|primers|sequence of reverse primer|
&nbsp;||    name| 806R|anything|primers|name of reverse primer|
primer_cutting|||||primers||arguments for primer cutting by cutadapt
&nbsp;|  overlap||10|1-length of primer|primers|minimum length of detected primer|
&nbsp;|  count||2|a positive integer|primers|maximum number of primers removed from each end|
&nbsp;|  filter_if_not_match||any|any or both|primers|reads are discarded if primer is not found on both or any end| any is the more strict setting; not used in single-end mode
&nbsp;|  perc_mismatch||0.2|0-1|primers|% mismatch between read and each primer|don't set this to 1
&nbsp;|  indels||"--no-indels"|"--no-indels" or ""|primers|whether indels in the primer sequence are allowed|
&nbsp;|  both_primers_in_read||false|false or true|primers|whether both primers are expected to be in the read| only used in single-end mode
sequencing_direction||| "unknown"|fwd_1, rvs_1 or unknown|primers| fwd_1: fwd primer in read 1; rvs_1: rvs primer in read 1; unknown: you don't know the sequencing direction or the direction is mixed |if you don't know the direction, dadasnake will try to re-orient using the primers
nextseq_novaseq||| false| true or false | primers| whether poly-G tails should be removed | set for Nextseq or Novaseq data
filtering |  | | | | dada | | settings for quality / length filtering; note on terminology: for paired sequencing fwd read refers to reads that had fwd primer or were declared as such (if no primer cutting was done); for single-end workflow, only the fwd setting is used, no matter the sequencing direction
&nbsp;|  trunc_length||||dada||length to truncate to (shorter reads are discarded)
&nbsp;||    fwd|0|a positive integer|dada|length after which fwd read is cut - shorter reads are discarded|0: no truncation by length; if you've cut the primers, this number refers to the length left after primer cutting
&nbsp;||    rvs|0|a positive integer|dada|length after which rvs read is cut - shorter reads are discarded|0: no truncation by length; ignored in single-ende mode; if you've cut the primers, this number refers to the length left after primer cutting
&nbsp;|  trunc_qual||||dada||length to truncate to (shorter reads are discarded)
&nbsp;| |    fwd|2|0-40|dada| fwd reads are cut before the first position with this quality| 
&nbsp;| |    rvs|2|0-40|dada| rvs reads are cut before the first position with this quality|ignored in single-ende mode
&nbsp;|  max_EE||||dada||filtering by maximum expected error after truncation: Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))
&nbsp;| |    fwd|2| a positive number |dada| After truncation, read pairs with higher than maxEE "expected errors" in fwd read will be discarded | use with trunc_length and/or truncQ; note that low truncQ or high trunc_length make it difficult to reach low maxEE values
&nbsp;||    rvs|2|a positive number|dada| After truncation, read pairs with higher than maxEE "expected errors" in rvs read will be discarded|ignored in single-ende mode; use with trunc_length and/or truncQ; note that low truncQ or high trunc_length make it difficult to reach low maxEE values
&nbsp;|  minLen||||dada||filtering by mimum length
&nbsp;||    fwd|20|a positive integer|dada|Remove reads with length less than minLen on fwd read. minLen is enforced after trimming and truncation.|use with truncQ
&nbsp;||    rvs|20|a positive integer|dada|Remove reads with length less than minLen on rvs read. minLen is enforced after trimming and truncation.| ignored in single-ende mode; use with truncQ
&nbsp;|  maxLen||||dada||filtering by maximum length
&nbsp;||    fwd| Inf|a positive integer or Inf|dada|Remove reads with length of fwd read greater than maxLen. maxLen is enforced before trimming and truncation.|
&nbsp;||    rvs| Inf|a positive integer or Inf|dada|Remove reads with length of rvs read greater than maxLen. maxLen is enforced before trimming and truncation.|ignored in single-ende mode
&nbsp;|  minQ||||dada||filtering by minimum quality after tuncation
&nbsp;||    fwd|0|0 or a positive number|dada|read pairs that contain a quality score lower than this in the fwd read after truncation will be discarded|use with trunc_length
&nbsp;||    rvs|0|0 or a positive number|dada|read pairs that contain a quality score lower than this in the rvs read after truncation will be discarded| ignored in single-ende mode; use with trunc_length
&nbsp;|  trim_left||||dada||
&nbsp;||    fwd|0|0 or a positive number|dada|this many bases will be cut from the 5' end of fwd reads|filtered reads will have length truncLen-trimLeft
&nbsp;||    rvs|0|0 or a positive number|dada|this many bases will be cut from the 5' end of rvs reads|filtered reads will have length truncLen-trimLeft
&nbsp;|  rm_phix||true|true or false|dada|remove phiX|useful with Illumina sequencing
downsampling|  | | | | dada | | 
&nbsp;|do||false|true or false|dada|set to true if you want to downsample before DADA2 ASV construction|
&nbsp;|number||50000|positive integer|dada|number of reads to keep per sample|
&nbsp;|min||true|true or false|dada|true to keep only samples with that many reads|samples with less reads are discarded
&nbsp;|use_total||false|true or false|dada|downsample to the fraction of a total number|useful for testing settings
&nbsp;|total||100000000|positive integer|dada|total number of reads to keep over all samples|used only with use_total
&nbsp;|seed||123|any positive integer|dada|seed for downsampling|keep constant in re-runs
error_seed|||100|any positive integer|dada|seed for error models|keep constant in re-runs
dada|||||dada||special DADA2 settings - default is good for Illumina
&nbsp;|  band_size||16|a positive integer|dada|Banding restricts the net cumulative number of insertion of one sequence relative to the other. | default is good for Illumina; set to 32 for 454 or PacBio
&nbsp;|  homopolymer_gap_penalty|| NULL|NULL or a negative integer|dada|The cost of gaps in homopolymer regions (>=3 repeated bases). Default is NULL, which causes homopolymer gaps to be treated as normal gaps.| default is good for Illumina; set to -1 for 454
&nbsp;|  pool|| false|true, false, "pseudo", or "within_run"|dada|Should DADA2 be run per sample (default) or in a pool, or should pseudo-pooling be done?| default is good for Illumina and much more efficient for large data sets; set to true for 454, pacbio and nanopore; set to pseudo for non-huge datasets, if you're interested in rare ASVs. You can also have within-run pools, but this setting is rarely useful.
&nbsp;|  omega_A|| 1e-40|number between 0 and 1|dada|Threshold to start new partition based on abundance in ASV finding.| default is good for Illumina; set lower for 454; according to the DADA2 authors, it's an underused feature - it can also kill your analysis
&nbsp;|  priors|| ""|"" or the absolute path to a fasta file with prior sequence data|dada|You can give DADA2 sequences to look out for in your dataset.| Don't change unless you know what you're doing.
&nbsp;|  omega_P|| 1e-4|number between 0 and 1|dada|Like omega_A, but for sequences matched by priors.| Only does anything, if you gave priors.
&nbsp;|  omega_C|| 1e-40|number between 0 and 1|dada|Threshold to start new partition based on quality in ASV finding.| Don't change unless you know what you're doing.
&nbsp;|  selfConsist|| false|true or false|dada|Should DADA2 do multiple rounds of ASV inference based on the normal error estimation?| Don't change unless you know what you're doing.
&nbsp;|  no_error_assumptions|| false|true or false|dada|If you've set selfConsist to true, you can make DADA2 not start from the normal error estimation.| Don't change unless you know what you're doing.
&nbsp;|  errorEstimationFunction|| loessErrfun|loessErrfun, PacBioErrfun or noqualErrfun, or loessErrfun_mod1 to 4|dada|The error estimation method within the DADA2 inference step.| default is good for Illumina; set to PacBioErrfun for pacbio and possibly to noqualErrfun if your hacking data without real quality values; ErnakovichLab models for Novaseq data are also available as e.g. loessErrfun_mod4
&nbsp;|  use_quals|| true|true or false|dada|DADA2 can be run without caring about quality.| Don't change unless you know what you're doing.
&nbsp;|  gapless|| true|true or false|dada|In the pre-screening, Kmers are employed to find gaps.| Don't change unless you know what you're doing - might help with 454 data and the like.
&nbsp;|  kdist_cutoff|| 0.42|a number between 0 and 1|dada|After the pre-screening, sequences of Kmers with this similarity are checked for actual matches.| Don't change unless you know what you're doing.
&nbsp;|  match|| 4|a number|dada|Score for match in Needleman-Wunsch-Alignment (the check for matching sequences).| Don't change unless you know what you're doing.
&nbsp;|  mismatch|| -5|a number|dada|Penaltiy for mismatch in Needleman-Wunsch-Alignment (the check for matching sequences).| Don't change unless you know what you're doing.
&nbsp;|  gap_penalty|| -8|a number|dada|Penaltiy for gaps in Needleman-Wunsch-Alignment (the check for matching sequences), unless the gaps are part of homopolymers - these are handled separately, see above.| Don't change unless you know what you're doing.
pair_merging|||||dada||settings for merging of read pairs
&nbsp;|  min_overlap||12|a positive integer|dada|The minimum length of the overlap required for merging the forward and reverse reads.|ignored in single-ende mode
&nbsp;|  max_mismatch||0|0 or a positive integer|dada|The maximum mismatches allowed in the overlap region.|ignored in single-ende mode
&nbsp;|  just_concatenate|| false|true or false|dada|whether reads should be concatenated rather than overlapped| ignored in single-ende mode; If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them.
&nbsp;|  trim_overhang|| true|true or false|dada|whether overhangs should be trimmed off after merging| ignored in single-ende mode; usually, overhangs should have been removed with the primer cutting step
chimeras|||||dada||settings for chimera removal
&nbsp;|  remove|| true|true or false|dada|whether chimeras should be removed|
&nbsp;|  method|| consensus|consensus, pooled or per-sample|dada|how chimeras are detected| consensus: samples are checked individually and sequences are removed by consensus; pooled: the samples are pooled and chimeras are inferred from pool; samples are checked individually and sequence counts of chimeras are set to 0 in individual samples
&nbsp;|  minFoldParentOverAbundance|| 2|a number > 1|dada|how overabundant do parents have to be to consider a read chimeric?| Should be higher for long amplicons (e.g. pacbio 3.5)
&nbsp;|  minParentAbundance|| 8|a number > 1|dada|how abundant do parents have to be to consider a read chimeric?| Don't change unless you know what you're doing.
&nbsp;|  allowOneOff|| false|true or false|dada|should sequences with a mismatch be flagged as potential chimera?| Don't change unless you know what you're doing.
&nbsp;|  minOneOffParentDistance|| 4|a number > 1|dada|if flagging sequences with one mismatch as potential one-off parents, how many mismatches are needed| Don't change unless you know what you're doing.
&nbsp;|  maxShift|| 16|a number|dada|maximum shift when aligning to potential parents| Don't change unless you know what you're doing.
taxonomy|||||taxonomy||settings for taxonomic annotation
&nbsp;|  dada||||taxonomy||settings for DADA2 implementation of bayesian classifier
&nbsp;||    do| false|true or false|taxonomy|whether DADA2 should be used for taxonomic annotation| the DADA2 implementation may work less well than the mothur classifier, and it may be slower
&nbsp;||    post_ITSx| false|true or false|taxonomy|whether the classifier should be run before or after ITSx| if you set this to true, you also have to set ITSx[do] to true; the DB isn't cut to a specific ITS region
&nbsp;||    run_on| list with ASV and cluster|list containing ASV and/or cluster|taxonomy|whether the classifier should be run on ASVs and or OTUs clustered from ASVs| 
&nbsp;||    db_path|"../DBs/DADA2"||taxonomy|directory where the database sits|change when setting up dadasnake on a new system
&nbsp;||    refFasta|"silva_nr99_v138_train_set.fa.gz"||taxonomy|training database name|
&nbsp;||    db_short_names|"silva_v138_nr99"||taxonomy|short name(s) to label database(s) in the output, separated by a whitespace; should be as many items as in ref_dbs_full|if your give less database names than databases, not all databases will be used
&nbsp;||    ref_dbs_full|""||taxonomy|full path and database file name(s) (without suffix), separated by a whitespace|if your give less database names than databases, not all databases will be used
&nbsp;||    minBoot|50|1-100|taxonomy|bootstrap value for classification|see DADA2 documentation for details
&nbsp;||    tryRC| false|false or true|taxonomy|if your reads are in the direction of the database (false), or reverse complement or you don't know (true)|true takes longer than false
&nbsp;||    seed|101|a positive integer|taxonomy|seed for DADA2 taxonomy classifier|keep constant in re-runs
&nbsp;||    look_for_species| false|true or false|taxonomy|whether you want to run a species-level annotation|species is an overkill for 16S data; if you set this, you need to have a specialised database (currently available for 16S silva 132)
&nbsp;||    spec_db|"../DBs/DADA2/silva_species_assignment_v138.fa.gz"||taxonomy|a DADA2-formatted species assignment database with path|change when setting up dadasnake on a new system
&nbsp;|  decipher||||taxonomy||settings for DECIPHER
&nbsp;||    do| false|true or false|taxonomy|whether DECIPHER should be used for taxonomic annotation| DECIPHER can work better than the mothur classifier, but it is slower and we don't have many databases for this software; you can run both DECIPHER and mothur (in parallel)
&nbsp;||    post_ITSx| false|true or false|taxonomy|whether DECIPHER should be run before or after ITSx| if you set this to true, you also have to set ITSx[do] to true; the DB isn't cut to a specific ITS region
&nbsp;||    run_on| list with ASV and cluster|list containing ASV and/or cluster|taxonomy|whether the classifier should be run on ASVs and or OTUs clustered from ASVs| 
&nbsp;||    db_path|"../DBs/decipher"||taxonomy|directory where the database sits|change when setting up dadasnake on a new system
&nbsp;||    tax_db|"SILVA_SSU_r138_2019.RData"||taxonomy|decipher database name|
&nbsp;||    db_short_names|"SILVA_138_SSU"||taxonomy|short name(s) to label database(s) in the output, separated by a whitespace; should be as many items as in ref_dbs_full|if your give less database names than databases, not all databases will be used
&nbsp;||    ref_dbs_full|""||taxonomy|full path and database file name(s) (without suffix), separated by a whitespace|if your give less database names than databases, not all databases will be used
&nbsp;||    threshold|60|1-100|taxonomy|threshold for classification|see DECIPHER documentation for details
&nbsp;||    strand| bottom|bottom, top or both|taxonomy|if your reads are in the direction of the database (top), reverse complement (bottom) or you don't know (both)|both takes roughly twice as long as the others
&nbsp;||    bootstraps|100|a positive integer|taxonomy|number of bootstraps|
&nbsp;||    seed|100|a positive integer|taxonomy|seed for DECIPHER run|keep constant in re-runs
&nbsp;||    look_for_species| false|true or false|taxonomy|whether you want to run a species-level annotation after DECIPHER|species is an overkill for 16S data; if you set this, you need to have a specialised database (currently available for 16S silva 132)
&nbsp;||    spec_db|"../DBs/DADA2/silva_species_assignment_v138.fa.gz"||taxonomy|a DADA2-formatted species assignment database with path|change when setting up dadasnake on a new system
&nbsp;|  mothur||||taxonomy||settings for Bayesian classifier (mothur implementation)
&nbsp;||    do| true|true or false|taxonomy|whether mothur's classify.seqs should be used for taxonomix annotation|we have more and more specific databases for mothur (and can make new ones), it's faster than DECIPHER, but potentially less correct; you can run both mothur and DECIPHER (in parallel)
&nbsp;||    post_ITSx| false|true or false|taxonomy|whether mothur's classify.seqs should be run before or after ITSx|if you set this to true, you also have to set ITSx[do] to true; use an ITSx-cut database if run afterwards
&nbsp;||    run_on| list with ASV and cluster|list containing ASV and/or cluster|taxonomy|whether the classifier should be run on ASVs and or OTUs clustered from ASVs| 
&nbsp;||    db_path|"../DBs/amplicon"||taxonomy|directory where the database sits|change when setting up dadasnake on a new system
&nbsp;||    tax_db|"SILVA_138_SSURef_NR99_prok.515F.806R"||taxonomy|the beginning of the filename of a mothur-formatted database|don't add .taxonomy or .fasta
&nbsp;||    db_short_names|"SILVA_138_SSU_NR99"||taxonomy|short name(s) to label database(s) in the output, separated by a whitespace; should be as many items as in ref_dbs_full|if your give less database names than databases, not all databases will be used
&nbsp;||    ref_dbs_full|""||taxonomy|full path and database file name(s) (without suffix), separated by a whitespace|if your give less database names than databases, not all databases will be used
&nbsp;||    cutoff|60|1-100|taxonomy|cut-off for classification|
blast|||||taxonomy||
&nbsp;|    do||true|true or false|taxonomy|whether blast should be run|
&nbsp;||    run_on| list with ASV and cluster|list containing ASV and/or cluster|taxonomy|whether blast should be run on ASVs and or OTUs clustered from ASVs| 
&nbsp;|    db_path||"../DBs/ncbi_16S_ribosomal_RNA"||taxonomy|path to blast database|
&nbsp;|    tax_db||16S_ribosomal_RNA||taxonomy|name (without suffix) of blast database|
&nbsp;|    e_val||0.01||taxonomy|e-value for blast|
&nbsp;|    tax2id||""|"tax2id table or "none"|taxonomy|whether taxonomic data is available in a tax2id table|this also assumes there is a taxdb file in the db_path; you don't need it, if you have a blast5 database
&nbsp;|    all||true||taxonomy|whether blastn should also be run on sequences that have been classified already|
&nbsp;|    run_basta||true|true or false|taxonomy|whether BASTA should be run on the BLASTn output|
&nbsp;|    basta_db||"../DBs/ncbi_taxonomy"||taxonomy|path to the NCBI-taxonomy database that is prepared when basta is installed|
&nbsp;|    basta_e_val||0.00001||taxonomy|e-value for hit selection|
&nbsp;|    basta_alen||100||taxonomy|minimum alignment length of hits|
&nbsp;|    basta_number||0|0 or a positive integer|taxonomy|maximum number of hits to use for classification|if set to 0 all hits will be considered
&nbsp;|    basta_min||3|a positive number|taxonomy|minimum number of hits a sequence must have to be assigned an LCA|needs to be smaller or equal to max_targets
&nbsp;|    basta_id||80|1-100|taxonomy|minimum identity of hit to be considered good|
&nbsp;|    basta_besthit||true|true or false|taxonomy|if set the final taxonomy will contain an additional column containing the taxonomy of the best (first) hit with defined taxonomy|
&nbsp;|    basta_perchits||99|an odd number greater than 50|taxonomy|percentage of hits that are used for LCA estimation|
ITSx|||||taxonomy||settings for ITSx
&nbsp;|  do|| false|true or false|taxonomy|whether ITSx should be run|only makes sense for analyses targetting an ITS region
&nbsp;|  min_regions||1|1-4|taxonomy|minimum number of detected regions|counting includes SSU, LSU and 5.8 next to the ITS regions
&nbsp;|  region|| ITS2|ITS1 or ITS2|taxonomy|which region to extract|
&nbsp;|  e_val||1.00E-05|0-1|taxonomy|e-value for ITS detection|
&nbsp;|  query_taxa||.|a letter|taxonomy|Profile set to use for the search|ITSx's -t option, see [manual](https://microbiology.se/publ/itsx_users_guide.pdf) for list
&nbsp;|  target_taxon||F|a letter|taxonomy|taxon output from ITSx to filter for|default is F for fungi
postclustering|||||dada||settings for clustering ASVs into OTUs (since 0.11)
&nbsp;|  do||true|true or false|dada|whether to do clustering, if no taxonomy is done|this is ignored if any of the taxonomy steps ask for clustered input
&nbsp;|  cutoff||0.97|a value between 0.5 and 1|dada|similarity cut-off|
&nbsp;|  method||vsearch|vsearch or deciperh|dada|clustering algorithm|
&nbsp;|  strand||plus|plus or both|dada|which strand to use for vsearch clustering|only used by vsearch, plus is faster and should be appropriate unless sequencing direction is unknown and can't be determined
final_table_filtering|||||postprocessing||settings for filtering the final ASV and/or OTU tables (before postprocessing, if postprocessing is done)
&nbsp;|do||true|true or false|postprocessing|whether a filtered version of the ASV/OTU table and sequences should be made and used for the post-processing steps|
&nbsp;|  keep_target_taxa||"."|"." or a regular expression for taxa to keep, e.g. "Bacteria"|postprocessing|pattern to look for in the taxstrings| done based on mothur and dada/DECIPHER result; "." means all are kept; all taxstrings are searched, if multiple classifiers were used - for clustered OTU tables, only the annotation of the OTUs is used, not the summary of ASV taxonomies
&nbsp;|target_min_length||0||postprocessing|minimal length sequence|
&nbsp;|target_max_length||Inf||postprocessing|maximum length of sequence|
postprocessing|||||postprocessing||settings for postprocessing
&nbsp;|  fungalTraits||||postprocessing||settings for fungalTraits
&nbsp;||    do|false|true or false|postprocessing|whether fungalTraits should be assigned|
&nbsp;||    db|"../DBs/functions/FungalTraits_1.2_ver_16Dec_2020_V.1.2.tsv"||postprocessing|path to fungalTraits DB|change when setting up dadasnake on a new system
&nbsp;||    classifier_db|mothur.SILVA_138_SSURef_NR99_cut||postprocessing|which classifier to use|can only be one
&nbsp;|  funguild||||postprocessing||settings for funguild
&nbsp;||    do|false|true or false|postprocessing|whether funguild should be run|
&nbsp;||    funguild_db|"../DBs/functions/funguild_db.json"||postprocessing|path to funguild DB|change when setting up dadasnake on a new system
&nbsp;||    classifier_db|mothur.SILVA_138_SSURef_NR99_cut||postprocessing|which classifier to use|can only be one
&nbsp;|  picrust2||||postprocessing||settings for PICRUSt2
&nbsp;||    do|true|true or false|postprocessing|whether PICRUSt2 should be run|
&nbsp;||    stratified|true|true or false|postprocessing|whether PICRUSt2 should return stratefied output|takes longer  
&nbsp;||    per_sequence_contrib|true|true or false|postprocessing|whether PICRUSt2 should run per_sequence_contrib routine|takes longer
&nbsp;||    skip_norm|false|true or false|postprocessing|whether PICRUSt2 should skip normalization of marker genes|
&nbsp;||    max_nsti|2|integer|postprocessing|PICRUSt2 max_nsti setting| see PICRUSt2 documentation for details  
&nbsp;||    do_nsti|true|true or false|postprocessing|whether PICRUSt2 should do NSTI| see PICRUSt2 documentation for details  
&nbsp;||    do_minpath|true|true or false|postprocessing|PICRUSt2 minpath setting| see PICRUSt2 documentation for details  
&nbsp;||    do_gapfill|true|true or false|postprocessing|PICRUSt2 gapfill setting| see PICRUSt2 documentation for details  
&nbsp;||    do_coverage|false|true or false|postprocessing|PICRUSt2 coverage setting| see PICRUSt2 documentation for details  
&nbsp;||    pathways|true|true or false|postprocessing|PICRUSt2 pathway setting| see PICRUSt2 documentation for details 
&nbsp;||    min_reads|1|integer|postprocessing|minimum number of reads per ASV to filter before PICRUSt2| setting this higher will remove rare ASVs from calculation (quicker and potentially less noisy)
&nbsp;||    min_samples|1|integer|postprocessing|minimum number of samples an ASV needs to be in before PICRUSt2| setting this higher will remove rare ASVs from calculation (quicker and potentially less noisy)
&nbsp;||    placement_tool|epa-ng|epa-ng or sepp|postprocessing|PICRUSt2 placement_tool setting| see PICRUSt2 documentation for details  
&nbsp;||    in_traits|EC,KO|comma-separated combination of COG, EC, KO, PFAM, TIGRFAM|postprocessing|PICRUSt2 in_traits setting| see PICRUSt2 documentation for details  
&nbsp;||    hsp_method|mp|mp, emp_prob, pic, scp, or subtree_average|postprocessing|PICRUSt2 hsp_method setting| see PICRUSt2 documentation for details 
&nbsp;||    edge_exponent|0.5|number|postprocessing|PICRUSt2 edge_exponent setting| see PICRUSt2 documentation for details  
&nbsp;||    min_align|0|number|postprocessing|PICRUSt2 min_align setting| see PICRUSt2 documentation for details
&nbsp;||    custom_trait_tables|''|string|postprocessing|PICRUSt2 custom_trait_tables setting| see PICRUSt2 documentation for details - not tested in dadasnake context yet
&nbsp;||    marker_gene_table|''|string|postprocessing|PICRUSt2 marker_gene_table setting| see PICRUSt2 documentation for details - not tested in dadasnake context yet
&nbsp;||    pathway_map|''|string|postprocessing|PICRUSt2 pathway_map setting| see PICRUSt2 documentation for details - not tested in dadasnake context yet
&nbsp;||    reaction_func|''|string|postprocessing|PICRUSt2 reaction_func setting| see PICRUSt2 documentation for details - not tested in dadasnake context yet
&nbsp;||    regroup_map|''|string|postprocessing|PICRUSt2 regroup_map setting| see PICRUSt2 documentation for details - not tested in dadasnake context yet
&nbsp;|  tax4fun2||||postprocessing||settings for tax4fun2 - deprecated !
&nbsp;||    do|false|true or false|postprocessing|whether tax4fun2 should be used|
&nbsp;||    db|"../DBs/functions/Tax4Fun2_ReferenceData_v2"||postprocessing|path to tax4fun2 DB|change when setting up dadasnake on a new system
&nbsp;||    database_mod|Ref99NR|Ref99NR or Ref100NR|postprocessing|which database to use|
&nbsp;||    normalize_by_copy_number|true|true or false|postprocessing|whether to normalize tax4fun2 results by copy number|normalization of pathway results is not possible
&nbsp;||    min_identity_to_reference|0.97|90 to 100 or 0.9 to 1.0|postprocessing|minimum similarity between ASV sequence and tax4fun DB|
&nbsp;||    user_data|false|true or false|postprocessing|whether user database should be used|
&nbsp;||    user_dir|"../DBs/Functions/GTDB_202_tax4fun2"||postprocessing|path to user database|
&nbsp;||    user_db|GTDB_fun||postprocessing|path to user database|
&nbsp;|  treeing|||true or false|postprocessing||
&nbsp;||    do|true||postprocessing|whether a phylogenetic tree should be made|
&nbsp;||    fasttreeMP|""||postprocessing|path to fasttreeMP executable|change when setting up dadasnake on a new system
&nbsp;|  rarefaction_curve||true|true or false|postprocessing|whether a rarefaction curve should be made|
sessionName |  |  | "" | "" or a single word | all | session name | only read, if you're not using the dadasnake wrappernormalMem |  |  | "" | "" or a number and letter | all | size of the RAM of one core of your normal copute nodes (e.g. 8G) | may be fixed during installation, only necessary for cluster submission 
bigMem |  |  | "" | "" or a number and letter | all | size of the RAM of one core of your high memory copute nodes (e.g. 30G) | may be fixed during installation, only necessary for cluster submission
bigCores |  |  | "" | "" or a number | all | maximum number of high memory copute nodes to use (e.g. 4) | 0 means all nodes have the same (normal) size may be fixed during installation, only necessary for cluster submission
sessionKind |  |  | "" | a string | all | automatically set by dadasnake wrapper | keep "" 
settingsLocked |  |  | false | a boolean or string | all | automatically set by dadasnake wrapper | it doesn't matter what you do


## What if something goes wrong?
If you gave dadasnake your email address and your system supports mailing (to that address), you will receive an email upon start and if the workflow encountered a problem or after the successful run. If there was a problem, you have to check the output and logs.
* Use the -d option of dadasnake or the --dryrun option of Snakemake before the run to check that your input files are where you want them and that you have permissions to write to your target directory. This will also do some checks on the configuration and samples table, so it discovers the majority of errors on a suitable combination of dataset and configuration.
* You can not make two runs of dadasnake write to the same output directory. If you start the second run while the first is still running, you will get an error either indicating that the directory can't be locked, or that the metadata is incomplete. If you've finished the first run already, the dadasnake will tell you that there's nothing to be done. Change the output directory in the config file to be unique for each run.
* A common reason for errors are misformatted inputs, e.g. the databases for the classification or the read files.
* dadasnake should catch most errors related to empty outputs. For example: the filtering is too stringent and no sequences are left; the primers you expected to find are not present; the sequences were truncated too short to be merged. Please report issues where this didn't happen.
* The best way to pinpoint those errors is to first check the .stderr file made by dadasnake (or the Snakemake output, if you run the workflow outside dadasnake). This will tell you which rule encountered the error, and, if you use the cluster submission, the job ID. You may have to search for the error a bit, because dadasnake will try to finish as much as possible of your run before dying. Hint: you can find errors by colour or by searching for "Error in rule".
* If you use the cluster submission, log files for every rule are written into the output directory and you can check the one with the job ID for additional information, otherwise the same information is written to the Snakemake output.
* The logs directory in the output directory contains log files for all steps that can produce comments. They are named with the step and then the name of the rule, so you can check the log file of the step that sent the error. Depending on the tool that sent the error, this will be easy to understand or cryptic. Don't hesitate to raise an issue in this repository if you get stuck.

## How to ...?
**I don't have primers on my reads, what do I do?**
Set `do_primers: false` in the configuration file, but make sure that orientation of the reads is the same.

**I did paired end sequencing, but my reads are too short to overlap**
You have two options: 
1) use only one read (usually the first) by setting `paired: false` in the config file and providing only the read you want to use in the samples table. This will run a single-end workflow. The makers of DADA2 would probably recommend this option in most cases.
2) use both reads, set a truncation length for filtering to make sure the sequences have the same lengths and use DADA2's option to "merge" reads without overlap e.g. 
```
filtering:
  trunc_length:
    fwd: 250
    rvs: 200
pair_merging:
  min_overlap: 0
  just_concatenate: true
```

**I need to set further parameters for job submission**
You can change the cluster configs and add the parameter, for example directly as part of the call field.

**I need to bind the jobs to the same node as the main job**
Yes, you can. If you use the submission-based wrapper, you can provide the flag for choosing a node as part of the SUBMIT_COMMAND variable in the VARIABLE_CONFIG file. Also, specify BIND_JOBS_TO_MAIN as true. You also need to set the variable that holds the node's name in your submission system as NODENAME_VAR. All jobs will then be submitted to the same node as the one that runs the main snakemake, if you include the flag for choosing a node as part of the call field in the cluster config. You can also specify that one, using -b. Example:
*VARIABLE_CONFIG* file:
```
...
SUBMIT_COMMAND	slurm --nodelist=
BIND_JOBS_TO_MAIN	true
NODENAME_VAR	SLURMD_NODENAME
SCHEDULER	slurm_simple
...
```

*slurm_simple.config*:
```
__default__:
  call: "sbatch --nodelist="
  mem_per_cpu: "--mem-per-cpu "
  partition: ""
  runtime: "-t"
  threads: "-c"
  stdout: "-o dadasnake.{rule}.{wildcards}.stdout"
```

*call*:
```
./dadasnake -c -b favorite_node -n TESTRUN config/config.test.yaml
```


**I have a very large dataset**
Great, if you have the computing power to match it, dadasnake will help you. It has successfully processed >27,000 samples in the same run. If you run out of memory in your run, set big_data to true in the config file and allow the use of multiple bigmem cores (we needed 360GB RAM for the 27,000 dataset). Disable highly memory intensive steps, such as treeing, chimera removal, plotting of rarefaction curves. If you didn't use the grouping by runs in your sample table, invent some runs of approx 100 samples each - these will be treated separately for some of the heavier DADA2 steps (error estimation).

**How do I restart a failed run?**
Depends on why it failed...
* If you ran into a time limit or similar, you can just run dadasnake on the same config with the -u option and then again with the -c option. This will make Snakemake pick up where it left off.
* For most other situations, it's probably best to fix what caused the error in your config file and delete the output directory to start from scratch. If you're going to be loosing a lot of run time to that, and you're quite certain the problem is only in the last attempted step, you can try to restart. Ask us, if in doubt.

**Can I restart from a certain step?**
If you're familiar with Snakemake, you can use it to force re-running the steps you need. It's not (yet) part of the dadasnake to do this more comfortably.

