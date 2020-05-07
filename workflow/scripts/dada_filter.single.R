log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(BiocParallel)
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}else{
    parallel <- FALSE
    register(SerialParam())
}
library(dada2)


# File parsing
fastq <- snakemake@input[[1]]
filt <- snakemake@output[[1]]
sampleName <- gsub("preprocessing/","run.",gsub("/","",gsub(".fastq","",fastq)))
filtpath <- gsub("/.+","",filt)

#output table with number of reads pre-filtering

print("filtering")

fastqFilter(fastq, 
            filt,
            truncLen=snakemake@config[['filtering']][['trunc_length']][['fwd']], 
            maxEE=snakemake@config[['filtering']][['max_EE']][['fwd']],
            maxLen=snakemake@config[['filtering']][['maxLen']][['fwd']],
            minLen=snakemake@config[['filtering']][['minLen']][['fwd']],
            minQ=snakemake@config[['filtering']][['minQ']][['fwd']],
            truncQ=snakemake@config[['filtering']][['trunc_qual']], 
            maxN=snakemake@config[['filtering']][['maxN']],
            rm.phix=as.logical(snakemake@config[['filtering']][['rm_phix']]),
            compress=TRUE, verbose=TRUE)

#output table with number of reads post-filtering
print("filtering done")