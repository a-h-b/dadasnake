log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

.libPaths(paste0(snakemake@config[['dada_lib']],"/R/library"))

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
#library(ggplot2)


# File parsing
filts <- unlist(snakemake@input)

sampleNames <- gsub("filtered/","run.",gsub("/","",gsub(".....fastq.gz","",filts)))

errfile <- snakemake@output[[1]]
dadafile <- snakemake@output[[3]]
errpath <- gsub("/.+","",errfile) 

print("filtering")

#if(!file_test("-d", errpath)) dir.create(errpath)
names(filts) <- sampleNames
set.seed(snakemake@config[['error_seed']])

print(paste0("learning error models ",snakemake@wildcards[['direction']]))
# Dereplication
dereps <- derepFastq(filts, verbose=TRUE)
# Name the derep-class objects by the sample names
names(dereps) <- sampleNames

errs <- learnErrors(filts, multithread=snakemake@threads)
saveRDS(errs,errfile)
pdf(snakemake@output[[2]],width=8,height=11,pointsize=7)
    plotErrors(errs, nominalQ=TRUE)
dev.off()

print("dada constructor")
dadas <- dada(dereps, err=errs, multithread=snakemake@threads)
saveRDS(dadas,dadafile)

print("done")
