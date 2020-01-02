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
#library(ggplot2)


# File parsing
filtFs <- sort(grep("fwd.fastq",unlist(snakemake@input),value=T))
filtRs <- sort(grep("rvs.fastq",unlist(snakemake@input),value=T))

if(length(filtFs) != length(filtRs)) stop("Forward and reverse files do not match.")

sampleNamesF <- gsub("filtered/","run.",gsub("/","",gsub(".fwd.fastq.gz","",filtFs)))
sampleNamesR <- gsub("filtered/","run.",gsub("/","",gsub(".rvs.fastq.gz","",filtRs)))

if(!all(sampleNamesF==sampleNamesR)) stop("Forward and reverse files do not match.")

errfileF <- snakemake@output[[1]]
errfileR <- snakemake@output[[2]]
dadafileF <- snakemake@output[[5]]
dadafileR <- snakemake@output[[6]]
errpath <- gsub("/.+","",errfileF) 

print("filtering")

#if(!file_test("-d", errpath)) dir.create(errpath)
names(filtFs) <- sampleNamesF
names(filtRs) <- sampleNamesR
set.seed(snakemake@config[['error_seed']])

print("learning error models F")
# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNamesF
names(derepRs) <- sampleNamesR

errFs <- learnErrors(filtFs, multithread=snakemake@threads)
saveRDS(errFs,errfileF)
pdf(snakemake@output[[3]],width=8,height=11,pointsize=7)
    plotErrors(errFs, nominalQ=TRUE)
dev.off()

print("learning error models R")
errRs <- learnErrors(filtRs, multithread=snakemake@threads)
saveRDS(errRs,errfileR)
pdf(snakemake@output[[4]],width=8,height=11,pointsize=7)
    plotErrors(errRs, nominalQ=TRUE)
dev.off()

print("dada constructor")
dadaFs <- dada(derepFs, err=errFs, multithread=snakemake@threads)
saveRDS(dadaFs,dadafileF)
dadaRs <- dada(derepRs, err=errRs, multithread=snakemake@threads)
saveRDS(dadaRs,dadafileR)

print("done")
