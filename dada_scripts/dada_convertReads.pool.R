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

filt <- unlist(snakemake@input[[2]])

filtNames <- sapply(filt,
                    function(x){
                      h <- unlist(strsplit(x,split="/"))
                      h <- h[length(h)-1:0]
                      paste(h,sep="/",collapse="/")
                    })

sampleNames <- gsub(".fastq.gz","",filtNames)

names(filt) <- sampleNames

mergefile <- snakemake@output[[1]]

print("merging")
err <- readRDS(errfile)

# Sample inference and merger of paired-end reads
derep <- derepFastq(filt)
if(snakemake@params[['pool']]=="pseudo"){
  print(paste0("make pseudo-pooled dada object"))
  dada <- dada(derep, err=err, multithread=snakemake@threads,
               BAND_SIZE=snakemake@config[['dada']][['band_size']],
               HOMOPOLYMER_GAP_PENALTY=snakemake@config[['dada']][['homopolymer_gap_penalty']],
               pool="pseudo")
}else{
  print(paste0("make pooled dada object"))
  dada <- dada(derep, err=err, multithread=snakemake@threads,
               BAND_SIZE=snakemake@config[['dada']][['band_size']],
               HOMOPOLYMER_GAP_PENALTY=snakemake@config[['dada']][['homopolymer_gap_penalty']],
               pool=TRUE)
  
}
saveRDS(dada,mergefile)
print("done")
