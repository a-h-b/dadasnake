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
The config file must be in .yaml format. The order within the yaml file does not matter, but the hierarchy has to be kept. Here are some explanations.

**top-level parameters** | **sub-parameters** | **subsub-parameters** | **default value** | **possible values** | **used in stage** | **explanation** | **comments / recommendations**
---|---|---|---|---|---|---|---
email |  | | "" | "" or a valid email address | all | email address for mail notification | keep empty if you don't want emails. Check spelling, it's not tested.
tmp_dir |  |  | "/work/$USER/tmp" | any path that you have permissions for writing to | all | directory for temporary, intermediate files that shouldn't be kept | keep this in your /work so you don't need to worry about removing its contents
raw_directory |  |  | "/work/$USER"| any one path where you might have your raw data | all | directory with all raw data | you will usually have this somewhere in a project folder
sample_table |  |  | "/work/$USER/samples.tsv" | any one location of your samples table | all | path to the samples table | you can keep this in your /work, because the dadasnake will copy it to your output directory
outputdir |  |  | "/data/project/metaamp/PLAYGROUND" | any path that you have permissions for writing to | all | directory where all the output will go | change this; a subdirectory of /work works best, but remember to move to a steady location afterwards; each output directory can hold the results of one completed pipeline only
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
paired||| true|true or false|primers and dada|do you want to use paired-end sequencing data?|if true, you have to give r1_file and r2_file in the samples table, if false only r1_file is read (if you want to use only R2 files from a paired-end sequencing run, put them in the r1_file column)
sequencing_direction||| "unknown"|fwd_1, rvs_1 or unknown|primers| fwd_1: fwd primer in read 1; rvs_1: rvs primer in read 1; unknown: you don't know the sequencing direction or (for paired-end sequencing) the direction is mixed |if you want to run single-end data and don't know the direction, you have to establish this first, because the dadasnake will not help you
primer_cutting|||||primers||arguments for primer cutting by cutadapt
&nbsp;|  overlap||10|1-length of primer|primers|minimum length of detected primer|
&nbsp;|  count||2|a positive integer|primers|maximum number of primers removed from each end|
&nbsp;|  filter_if_not_match||any|any or both|primers|reads are discarded if primer is not found on both or any end| any is the more strict setting; not used in single-end mode
&nbsp;|  perc_mismatch||0.2|0-1|primers|% mismatch between read and each primer|don't set this to 1
&nbsp;|  indels||"--no-indels"|"--no-indels" or ""|primers|whether indels in the primer sequence are allowed|
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
pair_merging|||||dada||settings for merging of read pairs
&nbsp;|  min_overlap||12|a positive integer|dada|The minimum length of the overlap required for merging the forward and reverse reads.|ignored in single-ende mode
&nbsp;|  max_mismatch||0|0 or a positive integer|dada|The maximum mismatches allowed in the overlap region.|ignored in single-ende mode
&nbsp;|  just_concatenate|| false|true or false|dada|whether reads should be concatenated rather than overlapped| ignored in single-ende mode; If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them.
&nbsp;|  trim_overhang|| true|true or false|dada|whether overhangs should be trimmed off after merging| ignored in single-ende mode; usually, overhangs should have been removed with the primer cutting step
chimeras|||||dada||settings for chimera removal
&nbsp;|  remove|| true|true or false|dada|whether chimeras should be removed|
&nbsp;|  method|| consensus|consensus, pooled or per-sample|dada|how chimeras are detected| consensus: samples are checked individually and sequences are removed by consensus; pooled: the samples are pooled and chimeras are inferred from pool; samples are checked individually and sequence counts of chimeras are set to 0 in individual samples
taxonomy|||||taxonomy||settings for taxonomic annotation
&nbsp;|  decipher||||taxonomy||settings for DECIPHER
&nbsp;||    do| true|true or false|taxonomy|whether DECIPHER should be used for taxonomic annotation| DECIPHER can work better than the mothur classifier, but it is slower and we don't have many databases for this software; you can run both DECIPHER and mothur (in parallel)
&nbsp;||    post_ITSx| false|true or false|taxonomy|whether DECIPHER should be run before or after ITSx| if you set this to true, you also have to set ITSx[do] to true; the DB isn't cut to a specific ITS region
&nbsp;||    db_path|"/data/project/metaamp/DBs/decipher"||taxonomy|directory where the database sits|don't change
&nbsp;||    tax_db|"SILVA_SSU_r132_March2018.RData"||taxonomy|decipher database name|
&nbsp;||    threshold|60|1-100|taxonomy|threshold for classification|see DECIPHER documentation for details
&nbsp;||    strand| bottom|bottom, top or both|taxonomy|if your reads are in the direction of the database (top), reverse complement (bottom) or you don't know (both)|both takes roughly twice as long as the others
&nbsp;||    bootstraps|100|a positive integer|taxonomy|number of bootstraps|
&nbsp;||    seed|100|a positive integer|taxonomy|seed for DECIPHER run|keep constant in re-runs
&nbsp;||    look_for_species| false|true or false|taxonomy|whether you want to run a species-level annotation after DECIPHER|species is an overkill for 16S data; if you set this, you need to have a specialised database (currently available for 16S silva 132)
&nbsp;||    spec_db|"/data/project/metaamp/DBs/decipher/silva_species_assignment_v132.fa.gz"||taxonomy|a DADA2-formatted species assignment database with path|
&nbsp;|  mothur||||taxonomy||settings for Bayesian classifier (mothur implementation)
&nbsp;||    do| true|true or false|taxonomy|whether mothur's classify.seqs should be used for taxonomix annotation|we have more and more specific databases for mothur (and can make new ones), it's faster than DECIPHER, but potentially less correct; you can run both mothur and DECIPHER (in parallel)
&nbsp;||    post_ITSx| false|true or false|taxonomy|whether mothur's classify.seqs should be run before or after ITSx|if you set this to true, you also have to set ITSx[do] to true; use an ITSx-cut database if run afterwards
&nbsp;||    db_path|"/data/project/metaamp/DBs/amplicon"||taxonomy|directory where the database sits|don't change
&nbsp;||    tax_db|"ifoh_515f.iroh_806r.silva_132"||taxonomy|the beginning of the filename of a mothur-formatted database|
&nbsp;||    cutoff|60|1-100|taxonomy|cut-off for classification|
blast|||||taxonomy||
&nbsp;|    do|true||true or false|taxonomy|whether blast should be run on all non-annotated sequences|
&nbsp;|    db_path|"/data/db/ncbi/blast/db/nt/2018-09-07"|||taxonomy|path to blast database|
&nbsp;|    tax_db|nt|||taxonomy|name (without suffix) of blast database|
&nbsp;|    e_val|0.01|||taxonomy|e-value for blast|
&nbsp;|    tax2id|"/data/project/metaamp/DBs/ncbi/tax2ID.sorted.tsv"||"tax2id table or "none"|taxonomy|whether taxonomic data is available in a tax2id table|this also assumes there is a taxdb file in the db_path
ITSx|||||taxonomy||settings for ITSx
&nbsp;|  do|| false|true or false|taxonomy|whether ITSx should be run|only makes sense for analyses targetting an ITS region
&nbsp;|  min_regions||1|1-4|taxonomy|minimum number of detected regions|counting includes SSU, LSU and 5.8 next to the ITS regions
&nbsp;|  region|| ITS2|ITS1 or ITS2|taxonomy|which region to extract|
&nbsp;|  e_val||1.00E-05|0-1|taxonomy|e-value for ITS detection|
hand_off|||||dada, taxonomy, postprocessing||settings deciding if additional formats should be given
&nbsp;|  biom||true|true or false|dada, taxonomy|whether a biome format output should be written|biome contains OTU table or OTU table and taxonomy (if taxonomy was run); biome table is never filtered
&nbsp;|  phyloseq||true|true or false|taxonomy, postprocessing|whether a phyloseq object should be returned|contains OTU table and taxonomy and tree (if each was run; if tree is run on pruned OTU table, phyloseq object contains filtered dataset)
final_table_filtering|||||postprocessing||settings for filtering the final OTU table (before postprocessing, if postprocessing is done)
&nbsp;|do||true|true or false|postprocessing|whether a filtered version of the OTU table and sequences should be made and used for the post-processing steps|
&nbsp;|  keep_target_taxa||"."|"." or a regular expression for taxa to keep, e.g. "Bacteria"|postprocessing|pattern to look for in the taxstrings| done based on mothur and DECIPHER result; "." means all are kept; both taxstrings are searched, if both classifiers were used
&nbsp;|  length_filter||||postprocessing||settings for length filter
&nbsp;||target_min_length|0||postprocessing|minimal length sequence|doesn't care for ITSx results
&nbsp;||target_max_length|Inf||postprocessing|maximum length of sequence|doesn't care for ITSx results
postprocessing|||||postprocessing||settings for postprocessinf
&nbsp;|  funguild||||postprocessing||settings for funguild
&nbsp;||    do|false|true or false|postprocessing|whether funguild should be run|
&nbsp;||    funguild_db|"/data/project/metaamp/DBs/amplicon/funguild_db.json"||postprocessing|path to funguild DB|don't change
&nbsp;||    classifier|mothur|mothur or decipher, depending on what was used|postprocessing|which classifier to use|can only be one
&nbsp;|  treeing|||true or false|postprocessing||
&nbsp;||    do|true||postprocessing|whether a phylogenetic tree should be made|
&nbsp;||    fasttreeMP|"/data/project/metaamp/TOOLS/FastTreeMP"||postprocessing|path to fasttreeMP executable|don't change
&nbsp;|  rarefaction_curve||true|true or false|postprocessing|whether a rarefaction curve should be made|


## The samples table
Every samples table needs sample names (under header library) and file names (just the names, the path should be in the config file under header r1_file and potentially r2_file). Since DADA2 estimates run-specific errors, it can be helpful to give run IDs (under header run). If several fastq files should end up in the same column of the OTU table, you can indicate this by giving these libraries the same sample name (under header sample). Libraries from different runs are combined in the final OTU table (example 1). Libraries from the same run are combined after primer-processing (example 2).
Example 1:
![overview](https://github.com/a-h-b/dadasnake/blob/master/documentation/samples_ex1.png)
Example 2:
![overview](https://github.com/a-h-b/dadasnake/blob/master/documentation/samples_ex2.png)
