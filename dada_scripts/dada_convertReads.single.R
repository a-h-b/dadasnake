log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(BiocParallel)
if (snakemake@threads > 1) {
    library("BiocParallel")
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}else{
    parallel <- FALSE
    register(SerialParam())
}
library(dada2)

# File parsing
errfile <- snakemake@input[[1]]

filt <- snakemake@input[[2]]

sampleName <- gsub(".+/","",gsub(".fastq.gz","",filt))

names(filt) <- sampleName

mergefile <- snakemake@output[[1]]

print("merging")
err <- readRDS(errfile)

# Sample inference and merger of paired-end reads
print(paste0("make dada object, ",sampleName))
derep <- derepFastq(filt)
dada <- dada(derep, err=err, multithread=snakemake@threads,
             BAND_SIZE=snakemake@config[['dada']][['band_size']],
             HOMOPOLYMER_GAP_PENALTY=snakemake@config[['dada']][['homopolymer_gap_penalty']])
saveRDS(dada,mergefile)

print("done")
