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
fastqF <- snakemake@input[[1]]
fastqR <- snakemake@input[[2]]
filtF <- snakemake@output[[1]]
filtR <- snakemake@output[[2]]
sampleName <- gsub("preprocessing/","run.",gsub("/","",gsub(".fwd.fastq","",fastqF)))
filtpath <- gsub("/.+","",filtF)

#output table with number of reads pre-filtering

print("filtering")

fastqPairedFilter(c(fastqF,fastqR), 
                    c(filtF,filtR),
                    truncLen=c(snakemake@config[['filtering']][['trunc_length']][['fwd']],
                               snakemake@config[['filtering']][['trunc_length']][['rvs']]), 
                    maxEE=c(snakemake@config[['filtering']][['max_EE']][['fwd']],
                               snakemake@config[['filtering']][['max_EE']][['rvs']]), 
                    maxLen=c(snakemake@config[['filtering']][['maxLen']][['fwd']],
                               snakemake@config[['filtering']][['maxLen']][['rvs']]),       
                    minLen=c(snakemake@config[['filtering']][['minLen']][['fwd']],
                               snakemake@config[['filtering']][['minLen']][['rvs']]),
                    minQ=c(snakemake@config[['filtering']][['minQ']][['fwd']],
                               snakemake@config[['filtering']][['minQ']][['rvs']]),
                    truncQ=snakemake@config[['filtering']][['trunc_qual']],
                    maxN=snakemake@config[['filtering']][['maxN']],
                    rm.phix=as.logical(snakemake@config[['filtering']][['rm_phix']]), 
                    compress=TRUE, verbose=TRUE)

#output table with number of reads post-filtering
print("filtering done")
