log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(BiocParallel)
#parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}else{
    parallel <- FALSE
    register(SerialParam())
}
#register(SerialParam())
library(dada2)


# File parsing
fastq <- snakemake@input[[1]]
#fastqR <- snakemake@input[[2]]
filt <- snakemake@output[[1]]
#filtR <- snakemake@output[[2]]
sampleName <- gsub("preprocessing/","run.",gsub("/","",gsub(".fastq","",fastq)))
filtpath <- gsub("/.+","",filt)

#output table with number of reads pre-filtering

print("filtering")

#if(!file_test("-d", filtpath)) dir.create(filtpath)
fastqFilter(fastq, 
            filt,
            truncLen=snakemake@config[['filtering']][['trunc_length']][['fwd']], 
            maxEE=snakemake@config[['filtering']][['max_EE']][['fwd']],
            maxLen=snakemake@config[['filtering']][['maxLen']][['fwd']],
            minLen=snakemake@config[['filtering']][['minLen']][['fwd']],
            minQ=snakemake@config[['filtering']][['minQ']][['fwd']],
            truncQ=snakemake@config[['filtering']][['trunc_qual']], 
            maxN=0, rm.phix=TRUE,
            compress=TRUE, verbose=TRUE)

#output table with number of reads post-filtering
print("filtering done")
