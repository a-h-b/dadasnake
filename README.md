![logo](https://github.com/a-h-b/dadasnake/blob/master/documentation/snake_all-trans.png)

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
At this point, you have all the scripts you need to run the workflow using snakemake, and you'd just need to get some data and databases (see point 7). If you want to use the **comfortable dadasnake wrapper**, follow the points 2-6. 

2) Adjust the file VARIABLE_CONFIG to your requirements:
* SNAKEMAKE_VIA_CONDA - set this to true, if you don't have snakemake in your path and want to install it via conda. Leave empty, if you don't need an additional snakemake.
* LOADING_MODULES - insert a bash command to load modules, if you need them to run conda. Leave empty, if you don't need to load a module.
* SUBMIT_COMMAND - insert the bash command you'll usually use to submit a job to your cluster to run on a single cpu for a few days. You only need this, if you want to have the snakemake top instance running in a submitted job. You also have the option to run it on the frontend via tmux. Leave empty, if you want to use the frontend version. You need to have [tmux](https://github.com/tmux/tmux/wiki) installed.
* SCHEDULER - insert the name of the scheduler you want to use (currently `slurm` or `uge`). This determines the cluster config given to snakemake, e.g. the cluster config file for slurm is config/slurm.config.yaml
* MAX_THREADS - set this to the maximum number of cores you want to be using in a run. If you don't set this, the default will be 50. Users can override this setting at runtime.

3) Decide how you want to run dadasnake:
Only do one of the two:
* if you want to submit the top snakemake call to the cluster:
```
cp auxiliary_files/dadasnake_allSubmit dadasnake
chmod 755 dadasnake
```
* if you want to keep snakemake on the frontend:
```
cp auxiliary_files/dadasnake_tmux dadasnake
chmod 755 dadasnake
```

4) **optional**: Install snakemake via conda:
If you want to use snakemake via conda (and you've set SNAKEMAKE_VIA_CONDA to true), install the environment, as [recommended by Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):
```
conda install -c conda-forge mamba
mamba create --prefix $PWD/conda/snakemake_env
conda activate $PWD/conda/snakemake_env
mamba install -c conda-forge -c bioconda snakemake
conda deactivate
```
Alternatively, if the above does not work, you can install an older snakemake version without mamba like so:
```
conda env create -f workflow/envs/snakemake_env.yml --prefix $PWD/conda/snakemake_env
```
Dadasnake will run with Snakemake version >= 5.9.1 and hasn't been tested with any previous versions.

5) Set permissions / PATH:
Dadasnake is meant to be used by multiple users. Set the permissions accordingly. I'd suggest to have read access for all files for the users plus execution rights for the dadasnake file and the .sh scripts in the subfolder dada_scripts and the conda environment. Add the dadasnake directory to your path. It can also be useful to make the VARIABLE_CONFIG file not-writable, because you will always need it. The same goes for config.default.yaml once you've set the paths to the databases you want to use (see below).

6) Test run:
The test run does not need any databases. You should be able to start it by running 
```
./dadasnake -l -n "TESTRUN" -r config/config.test.yaml
```
If all goes well, dadasnake will run in the current session, load the conda environment, and make and fill a directory called testoutput. The first step (getting the conda environment) will take several minutes. A completed run contains a file "workflow.done". Don't worry if you see a few warnings from `mv`, such as `mv: cannot stat ‘slurm*’: No such file or directory`. 
If you don't want to see dadasnake's guts at this point, you can also run this with the -c or -f settings to submit to your cluster or start a tmux session (see How to run dadasnake below). 

The first run will install the conda environment containing DADA2 and the other programs that will be used by all users. I'd strongly suggest to **remove one line from the activation script** ( conda/XXXXXXXX/etc/conda/activate.d/activate-r-base.sh ) after the installation, namely the one reading: `R CMD javareconf > /dev/null 2>&1 || true`, because you don't need this line later and if two users run this at the same time it can cause trouble.

7) Databases:
The dadasnake does not supply databases. I'd suggest to use the [SILVA database](https://www.arb-silva.de/no_cache/download/archive/current/Exports/) for 16S data and [UNITE](https://doi.org//10.15156/BIO/786336) for ITS. 
* dadasnake uses [mothur](https://www.mothur.org/) to do the classification, as it's faster and likely more accurate than the DADA2 option. You need to format the database like for mothur ([see here](https://www.mothur.org/wiki/Taxonomy_outline)). 
* In addition to mothur, dadasnake implements [DECIPHER](http://www2.decipher.codes/Documentation.html). You can find decipher [databases](http://www2.decipher.codes/Downloads.html) on the decipher website or build them yourself. 
* You can also use dadasnake to blast and to annotate fungal taxonomy with guilds via funguild, if you have suitable databases. 
**You need to set the path to the databases of your choice in the config file.** By default, dadasnake looks for databases in the directory above where it was called. It makes sense to change this for your system in the config.default.yaml file upon installation, if all users access databases in the same place.

8) R-package phyloseq:
While DADA2 and other useful R-packages are part of the conda-environment, phyloseq does not like being installed via conda right now. If you want a phyloseq hand-off, install phyloseq into the common conda environment after the testrun. First, check which environment was created: 
```
ls ./conda
```
This will show a .yaml file and a directory with the same name. Use this name in activating the environment:
```
conda activate ./conda/XXXXXXXX
```
In the conda environment, run R:
```
R
```
Within R, set the path to the R-library in the conda environment and install phyloseq. This depends on regular building utilities of your operating system, e.g. gcc. Choose a mirror from the list when prompted. Don't update packages in the end (choose n):
```
.libPaths(paste0(Sys.getenv("CONDA_PREFIX"),"/lib/R/library"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")
```
Leave R without saving the workspace (choose n when prompted):
```
quit()
```
Deactivate the environment:
```
conda deactivate
```

9) Fasttree:
dadasnake comes with fasttree for treeing, but if you have a decent number of sequences, it is likely to be relatively slow. If you have fasttreeMP, you can give the path to it in the config file.

##  

![overview](https://github.com/a-h-b/dadasnake/blob/master/documentation/pipeline.png)

## How to run dadasnake
To run the dadasnake, you need a config file and a sample table, plus data: 
* The config file (in yaml format) is read by Snakemake to determine the inputs, steps, arguments and outputs. 
* The sample table (tab-separated text) gives sample names and file names in the simplest case, with column headers named library and r1_file (and r2_file for paired-end data sets) - its path has to be mentioned in the config file. You can add columns labeled `run` and `sample` to indicate libraries that should be combined into one final column and different sequencing runs. 
* All raw data (usually fastq files) need to be in one directory (which has to be given in the config file). 
* It is possible (and the best way to do this) to have a config file per run, which defines all settings that differ from the default config file.

### Using the dadasnake wrapper
As shown in the installation description above, dadasnake can be run in a single step, by calling dadasnake. Since most of the configuration is done via the config file, the options are very limited. You can either:
* -c run (submit to a cluster) dadasnake and make a report (-r), or
* -l run (in the current terminal) dadasnake and make a report (-r), or
* -f run (in a tmux session on the frontend) dadasnake *only available in the tmux installation* and make a report (-r), or
* just make a report (-r), or 
* run a dryrun (-d), or 
* unlock a working directory, if a run was killed (-u). 
It is strongly recommended to **first run a dryrun on a new configuration**, which will tell you within a few seconds and without submission to a cluster whether your chosen steps work together, the input files are where you want them, and your sample file is formatted correctly. In all cases you need the config file as the last argument. 
```
dadasnake -d -r config.yaml
```
You can add a name for your main job (-n NAME), e.g.:
```
dadasnake -c -n RUNNAME -r config.yaml
```
You can also set the number of cpus to maximally run at the same time with -t. The defaults (1 for local/frontend runs and 50 for clusters) are reasonable for many settings and if you don't know what this means, you probably don't have to worry. But you may want to increase the numbers for larger datasets or bigger infrastructure, or decrease the numbers to match your environment's constraints.
Depending on your dataset and settings and your cluster's scheduler, the workflow will take a few minutes to days to finish. 

If you use the tmux version, you can see the tmux process running by typing `tmux ls`. You can also see the progress by checking the stdandard error file `tail RUNNAME_XXXXXXXXXX.stderr`.


### Running snakemake manually
Once raw data, config file and sample file are present, the workflow can be started from the dadasnake directory by the snakemake command:
```
snakemake -s Snakefile --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda
```
If you're using a computing cluster, add your cluster's submission command and the number of jobs you want to maximally run at the same time, e.g.:
```
snakemake -j 50 -s Snakefile --cluster "qsub -l h_rt={params.runtime},h_vmem={params.mem} -pe smp {threads} -cwd" --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda
```
This will submit most steps as their own job to your cluster's queue. The same can be achieved with a [cluster configuration](https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html#cluster-execution):
```
snakemake -j 50 -s Snakefile --cluster-config PATH/TO/SCHEDULER.config.yaml --cluster "{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu} {cluster.threads}{threads} {cluster.partition}" --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda
```
If you want to share the conda installation with colleagues, use the `--conda-prefix` argument of Snakemake
```
snakemake -j 50 -s Snakefile --cluster-config PATH/TO/SCHEDULER.config.yaml --cluster "{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu} {cluster.threads}{threads} {cluster.partition}" --use-conda --conda-prefix /PATH/TO/YOUR/COMMON/CONDA/DIRECTORY
```
Depending on your dataset and settings, and your cluster's queue, the workflow will take a few minutes to days to finish.

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
You can control the settings for each step in a config file.

![steps](https://github.com/a-h-b/dadasnake/blob/master/documentation/steps.png)

## The configuration
The config file must be in .yaml format. The order within the yaml file does not matter, but the hierarchy has to be kept. Here are some explanations.

**top-level parameters** | **sub-parameters** | **subsub-parameters** | **default value** | **possible values** | **used in stage** | **explanation** | **comments / recommendations**
---|---|---|---|---|---|---|---
email |  | | "" | "" or a valid email address | all | email address for mail notification | keep empty if you don't want emails. Check spelling, it's not tested.
tmp_dir |  |  | "/work/$USER/tmp" | any path that you have permissions for writing to | all | directory for temporary, intermediate files that shouldn't be kept | keep this in your /work so you don't need to worry about removing its contents
raw_directory |  |  | "/work/$USER"| any one path where you might have your raw data | all | directory with all raw data | you will usually have this somewhere in a project folder
sample_table |  |  | "/work/$USER/samples.tsv" | any one location of your samples table | all | path to the samples table | you can keep this in your /work, because the dadasnake will copy it to your output directory
outputdir |  |  | "dadasnake_output" | any path that you have permissions for writing to | all | directory where all the output will go | change this; a scratch-type place works best (e.g. subdirectory of /work/$USER), but remember to move to a steady location afterwards; each output directory can hold the results of one completed pipeline run only
do_primers |  |  | true | true or false | all | should primers be cut? | 
do_dada |  |  | true | true or false | all | should DADA2 be run? | 
do_taxonomy |  |  | true | true or false | all | should taxonomic classification be done? |
do_postprocessing |  |  | true | true or false | all | should some more steps be done (e.g. functional annotation) | 
primers |  |  |  |  | primers |  | information on primers
  &nbsp;| fwd |  |  |  | primers |  | information on forward primer
  &nbsp;|  | sequence | GTGYCAGCMGCCGCGGTAA | any sequence of IUPAC DNA code | primers | sequence of forward primer |
 &nbsp;|  | name | 515F | anything | primers | name of forward primer | for your reference only
&nbsp;| rvs |  |  |  | primers |  | information on reverse primer
&nbsp;||    sequence| GGACTACNVGGGTWTCTAAT|any sequence of IUPAC DNA code|primers|sequence of reverse primer|
&nbsp;||    name| 806R|anything|primers|name of reverse primer|
paired||| true|true or false|primers and dada|do you want to use paired-end sequencing data?|if true, you have to give r1_file and r2_file in the samples table, if false only r1_file is read (if you want to use only R2 files from a paired-end sequencing run, put their name in the r1_file column)
sequencing_direction||| "unknown"|fwd_1, rvs_1 or unknown|primers| fwd_1: fwd primer in read 1; rvs_1: rvs primer in read 1; unknown: you don't know the sequencing direction or the direction is mixed |if you want to run single-end data and don't know the direction, dadasnake will re-orient the primers
primer_cutting|||||primers||arguments for primer cutting by cutadapt
&nbsp;|  overlap||10|1-length of primer|primers|minimum length of detected primer|
&nbsp;|  count||2|a positive integer|primers|maximum number of primers removed from each end|
&nbsp;|  filter_if_not_match||any|any or both|primers|reads are discarded if primer is not found on both or any end| any is the more strict setting; not used in single-end mode
&nbsp;|  perc_mismatch||0.2|0-1|primers|% mismatch between read and each primer|don't set this to 1
&nbsp;|  indels||"--no-indels"|"--no-indels" or ""|primers|whether indels in the primer sequence are allowed|
&nbsp;|  both_primers_in_read||false|false or true|primers|whether both primers are expected to be in the read| only used in single-end mode
filtering |  | | | | dada | | settings for quality / length filtering; note on terminology: for paired sequencing fwd read refers to reads that had fwd primer or were declared as such (if no primer cutting was done); for single-end workflow, only the fwd setting is used, no matter the sequencing direction
&nbsp;|  trunc_length||||dada||length to truncate to (shorter reads are discarded)
&nbsp;||    fwd|0|a positive integer|dada|length after which fwd read is cut - shorter reads are discarded|0: no truncation by length; if you've cut the primers, this number refers to the length left after primer cutting
&nbsp;||    rvs|0|a positive integer|dada|length after which rvs read is cut - shorter reads are discarded|0: no truncation by length; ignored in single-ende mode; if you've cut the primers, this number refers to the length left after primer cutting
&nbsp;|  trunc_qual||13|0-40|dada|reads are cut before the first position with this quality|
&nbsp;|  max_EE||||dada||filtering by maximum expected error after truncation: Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))
&nbsp;| |    fwd | 2 | a positive number | dada | After truncation, read pairs with higher than maxEE "expected errors" in fwd read will be discarded | use with trunc_length and/or truncQ; note that low truncQ or high trunc_length make it difficult to reach low maxEE values
&nbsp;||    rvs|2|a positive number|dada| After truncation, read pairs with higher than maxEE "expected errors" in rvs read will be discarded|ignored in single-ende mode; use with trunc_length and/or truncQ; note that low truncQ or high trunc_length make it difficult to reach low maxEE values
&nbsp;|  minLen||||dada||filtering by mimum length
&nbsp;||    fwd|200|a positive integer|dada|Remove reads with length less than minLen on fwd read. minLen is enforced after trimming and truncation.|use with truncQ
&nbsp;||    rvs|100|a positive integer|dada|Remove reads with length less than minLen on rvs read. minLen is enforced after trimming and truncation.| ignored in single-ende mode; use with truncQ
&nbsp;|  maxLen||||dada||filtering by maximum length
&nbsp;||    fwd| Inf|a positive integer or Inf|dada|Remove reads with length of fwd read greater than maxLen. maxLen is enforced before trimming and truncation.|
&nbsp;||    rvs| Inf|a positive integer or Inf|dada|Remove reads with length of rvs read greater than maxLen. maxLen is enforced before trimming and truncation.|ignored in single-ende mode
&nbsp;|  minQ||||dada||filtering by minimum quality after tuncation
&nbsp;||    fwd|0|0 or a positive number|dada|read pairs that contain a quality score lower than this in the fwd read after truncation will be discarded|use with trunc_length
&nbsp;||    rvs|0|0 or a positive number|dada|read pairs that contain a quality score lower than this in the rvs read after truncation will be discarded| ignored in single-ende mode; use with trunc_length
error_seed|||100|any positive integer|dada|seed for error models|keep constant in re-runs
dada|||||dada||special DADA2 settings - default is good for Illumina
&nbsp;|  band_size||16|a positive integer|dada|Banding restricts the net cumulative number of insertion of one sequence relative to the other. | default is good for Illumina; set to 32 for 454 or PacBio
&nbsp;|  homopolymer_gap_penalty|| NULL|NULL or a negative integer|dada|The cost of gaps in homopolymer regions (>=3 repeated bases). Default is NULL, which causes homopolymer gaps to be treated as normal gaps.| default is good for Illumina; set to -1 for 454
&nbsp;|  pool|| false|true, false or "pseudo"|dada|Should DADA2 be run per sample (default) or in a pool, or should pseudo-pooling be done?| default is good for Illumina and much more efficient for large data sets; set to true for 454, pacbio and nanopore; set to pseudo for non-huge datasets, if you're interested in rare ASVs.
&nbsp;|  omega_A|| 1e-40|number between 0 and 1|dada|Threshold to start new partition based on abundance in ASV finding.| default is good for Illumina; set lower for 454; according to the DADA2 authors, it's an underused feature - it can also kill your analysis
&nbsp;|  priors|| ""|"" or the absolute path to a fasta file with prior sequence data|dada|You can give DADA2 sequences to look out for in your dataset.| Don't change unless you know what you're doing.
&nbsp;|  omega_P|| 1e-4|number between 0 and 1|dada|Like omega_A, but for sequences matched by priors.| Only does anything, if you gave priors.
&nbsp;|  omega_C|| 1e-40|number between 0 and 1|dada|Threshold to start new partition based on quality in ASV finding.| Don't change unless you know what you're doing.
&nbsp;|  selfConsist|| false|true or false|dada|Should DADA2 do multiple rounds of ASV inference based on the normal error estimation?| Don't change unless you know what you're doing.
&nbsp;|  no_error_assumptions|| false|true or false|dada|If you've set selfConsist to true, you can make DADA2 not start from the normal error estimation.| Don't change unless you know what you're doing.
&nbsp;|  errorEstimationFunction|| loessErrfun|loessErrfun, PacBioErrfun or noqualErrfun|dada|The error estimation method within the DADA2 inference step.| default is good for Illumina; set to PacBioErrfun for pacbio and possibly to noqualErrfun if your hacking data without real quality values
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
&nbsp;|  decipher||||taxonomy||settings for DECIPHER
&nbsp;||    do| false|true or false|taxonomy|whether DECIPHER should be used for taxonomic annotation| DECIPHER can work better than the mothur classifier, but it is slower and we don't have many databases for this software; you can run both DECIPHER and mothur (in parallel)
&nbsp;||    post_ITSx| false|true or false|taxonomy|whether DECIPHER should be run before or after ITSx| if you set this to true, you also have to set ITSx[do] to true; the DB isn't cut to a specific ITS region
&nbsp;||    db_path|"../DBs/decipher"||taxonomy|directory where the database sits|change when setting up dadasnake on a new system
&nbsp;||    tax_db|"SILVA_SSU_r132_March2018.RData"||taxonomy|decipher database name|
&nbsp;||    threshold|60|1-100|taxonomy|threshold for classification|see DECIPHER documentation for details
&nbsp;||    strand| bottom|bottom, top or both|taxonomy|if your reads are in the direction of the database (top), reverse complement (bottom) or you don't know (both)|both takes roughly twice as long as the others
&nbsp;||    bootstraps|100|a positive integer|taxonomy|number of bootstraps|
&nbsp;||    seed|100|a positive integer|taxonomy|seed for DECIPHER run|keep constant in re-runs
&nbsp;||    look_for_species| false|true or false|taxonomy|whether you want to run a species-level annotation after DECIPHER|species is an overkill for 16S data; if you set this, you need to have a specialised database (currently available for 16S silva 132)
&nbsp;||    spec_db|"../DBs/decipher/silva_species_assignment_v132.fa.gz"||taxonomy|a DADA2-formatted species assignment database with path|change when setting up dadasnake on a new system
&nbsp;|  mothur||||taxonomy||settings for Bayesian classifier (mothur implementation)
&nbsp;||    do| true|true or false|taxonomy|whether mothur's classify.seqs should be used for taxonomix annotation|we have more and more specific databases for mothur (and can make new ones), it's faster than DECIPHER, but potentially less correct; you can run both mothur and DECIPHER (in parallel)
&nbsp;||    post_ITSx| false|true or false|taxonomy|whether mothur's classify.seqs should be run before or after ITSx|if you set this to true, you also have to set ITSx[do] to true; use an ITSx-cut database if run afterwards
&nbsp;||    db_path|"../DBs/amplicon"||taxonomy|directory where the database sits|change when setting up dadasnake on a new system
&nbsp;||    tax_db|"ifoh_515f.iroh_806r.silva_132"||taxonomy|the beginning of the filename of a mothur-formatted database|don't add .taxonomy or .fasta
&nbsp;||    cutoff|60|1-100|taxonomy|cut-off for classification|
blast|||||taxonomy||
&nbsp;|    do|false||true or false|taxonomy|whether blast should be run|
&nbsp;|    db_path|"../DBs/nt/2018-09-07"|||taxonomy|path to blast database|
&nbsp;|    tax_db|nt|||taxonomy|name (without suffix) of blast database|
&nbsp;|    e_val|0.01|||taxonomy|e-value for blast|
&nbsp;|    tax2id|""||"tax2id table or "none"|taxonomy|whether taxonomic data is available in a tax2id table|this also assumes there is a taxdb file in the db_path; you don't need it, if you have a blast5 database
&nbsp;|    all|false|||taxonomy|whether blastn should also be run on sequences that have been classified already|default means blast for run only on non-assigned sequences
ITSx|||||taxonomy||settings for ITSx
&nbsp;|  do|| false|true or false|taxonomy|whether ITSx should be run|only makes sense for analyses targetting an ITS region
&nbsp;|  min_regions||1|1-4|taxonomy|minimum number of detected regions|counting includes SSU, LSU and 5.8 next to the ITS regions
&nbsp;|  region|| ITS2|ITS1 or ITS2|taxonomy|which region to extract|
&nbsp;|  e_val||1.00E-05|0-1|taxonomy|e-value for ITS detection|
hand_off|||||dada, taxonomy, postprocessing||settings deciding if additional formats should be given
&nbsp;|  biom||true|true or false|dada, taxonomy|whether a biome format output should be written|biome contains OTU table or OTU table and taxonomy (if taxonomy was run); biome table is never filtered
&nbsp;|  phyloseq||false|true or false|taxonomy, postprocessing|whether a phyloseq object should be returned|contains OTU table and taxonomy and tree (if each was run; if tree is run on pruned OTU table, phyloseq object contains filtered dataset)
final_table_filtering|||||postprocessing||settings for filtering the final OTU table (before postprocessing, if postprocessing is done)
&nbsp;|do||true|true or false|postprocessing|whether a filtered version of the OTU table and sequences should be made and used for the post-processing steps|
&nbsp;|  keep_target_taxa||"."|"." or a regular expression for taxa to keep, e.g. "Bacteria"|postprocessing|pattern to look for in the taxstrings| done based on mothur and DECIPHER result; "." means all are kept; both taxstrings are searched, if both classifiers were used
&nbsp;|target_min_length||0||postprocessing|minimal length sequence|doesn't care for ITSx results
&nbsp;|target_max_length||Inf||postprocessing|maximum length of sequence|doesn't care for ITSx results
postprocessing|||||postprocessing||settings for postprocessinf
&nbsp;|  funguild||||postprocessing||settings for funguild
&nbsp;||    do|false|true or false|postprocessing|whether funguild should be run|
&nbsp;||    funguild_db|"../DBs/amplicon/funguild_db.json"||postprocessing|path to funguild DB|change when setting up dadasnake on a new system
&nbsp;||    classifier|mothur|mothur or decipher, depending on what was used|postprocessing|which classifier to use|can only be one
&nbsp;|  treeing|||true or false|postprocessing||
&nbsp;||    do|true||postprocessing|whether a phylogenetic tree should be made|
&nbsp;||    fasttreeMP|"FastTreeMP"||postprocessing|path to fasttreeMP executable|change when setting up dadasnake on a new system
&nbsp;|  rarefaction_curve||true|true or false|postprocessing|whether a rarefaction curve should be made|


## The samples table
Every samples table needs sample names (under header library) and file names (just the names, the path should be in the config file under header r1_file and potentially r2_file). Since DADA2 estimates run-specific errors, it can be helpful to give run IDs (under header run). If several fastq files should end up in the same column of the OTU table, you can indicate this by giving these libraries the same sample name (under header sample). Libraries from different runs are combined in the final OTU table (example 1). Libraries from the same run are combined after primer-processing (example 2).
Example 1:
![overview](https://github.com/a-h-b/dadasnake/blob/master/documentation/samples_ex1.png)
Example 2:
![overview](https://github.com/a-h-b/dadasnake/blob/master/documentation/samples_ex2.png)

## What if something goes wrong?
If you gave dadasnake your email address and your system supports mailing (to that address), you will receive an email upon start and if the workflow encountered a problem or after the successful run. If there was a problem, you have to check the output and logs.
* Use the -d option of dadasnake or the --dryrun option of Snakemake before the run to check that your input files are where you want them and that you have permissions to write to your target directory. This will also do some checks on the configuration and samples table, so it discovers the majority of errors on a suitable combination of dataset and configuration.
* You can not make two runs of dadasnake write to the same output directory. If you start the second run while the first is still running, you will get an error either indicating that the directory can't be locked, or that the metadata is incomplete. If you've finished the first run already, the dadasnake will tell you that there's nothing to be done. Change the output directory in the config file to be unique for each run.
* Reasons for dadasnake to encounter problems during the runs usually have something to do with empty outputs. For example: the filtering is too stringent and no sequences are left; the primers you expected to find are not present; the sequences were truncated too short to be merged. 
* The second reason are misformatted inputs, e.g. the databases for the classification or the read files.
* The best way to pinpoint those errors is to first check the .stderr file made by dadasnake (or the Snakemake output, if you run the workflow outside dadasnake). This will tell you which rule encountered the error, and, if you use the cluster submission, the job ID.
* If you use the cluster submission, log files for every rule are written into the output directory and you can check the one with the job ID for additional information, other wise the same information is written to the Snakemake output.
* The logs directory in the output directory contains log files for all steps that can produce comments. They are named with the step and then the name of the rule, so you can check the log file of the step that sent the error. Depending on the tool that sent the error, this will be easy to understand or cryptic. Don't hesitate to raise an issue at this repository if you get stuck.

## How to ...?
**I don't have primers on my reads, what do I do?**
Set `do_primers: false` in the configuration file, but make sure that their orientation is the same.

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

**How do I restart a failed run?**
Depends on why it failed...
* If you ran into a time limit or similar, you can just run dadasnake on the same config with the -u option and then again with the -c option. This will make Snakemake pick up where it left off.
* For most other situations, it's probably best to fix what caused the error in your config file and delete the output directory to start from scratch. If you're going to be loosing a lot of run time to that, and you're quite certain the problem is only in the last attempted step, you can try to restart. Ask us, if in doubt.

**Can I restart from a certain step?**
If you're familiar with Snakemake, you can use it to force re-running the steps you need. It's not (yet) part of the dadasnake to do this more comfortably.

