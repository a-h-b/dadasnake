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
library(Biostrings)

# File parsing
errfile <- snakemake@input[[1]]

filt <- unlist(snakemake@input[-1])

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

if(snakemake@config[['dada']][['priors']]!=""){
  priorFasta <- readDNAStringSet(snakemake@config[['dada']][['priors']])
  priors <- as.character(priorFasta)
}else{
  priors=""
}

if(as.logical(snakemake@config[['dada']][['no_error_assumptions']])){
  err <- NULL
}else{
  err <- readRDS(errfile)
}

# Sample inference and merger of paired-end reads
derep <- derepFastq(filt)
if(snakemake@params[['pooling']]=="pseudo"){
  print(paste0("make pseudo-pooled dada object"))
  dada <- dada(derep, err=err, multithread=snakemake@threads,
               BAND_SIZE=snakemake@config[['dada']][['band_size']],
               HOMOPOLYMER_GAP_PENALTY=snakemake@config[['dada']][['homopolymer_gap_penalty']],
               OMEGA_A=as.numeric(snakemake@config[['dada']][['omega_A']]),
               OMEGA_P=as.numeric(snakemake@config[['dada']][['omega_P']]),
               OMEGA_C=as.numeric(snakemake@config[['dada']][['omega_C']]),
               GAPLESS=as.logical(snakemake@config[['dada']][['gapless']]),
               KDIST_CUTOFF=as.numeric(snakemake@config[['dada']][['kdist_cutoff']]),
               MATCH=as.numeric(snakemake@config[['dada']][['match']]),
               MISMATCH=as.numeric(snakemake@config[['dada']][['mismatch']]),
               GAP_PENALTY=as.numeric(snakemake@config[['dada']][['gap_penalty']]),
               pool="pseudo",
               priors=priors,
               selfConsist=as.logical(snakemake@config[['dada']][['selfConsist']]),
               errorEstimationFunction=match.fun(snakemake@config[['dada']][['errorEstimationFunction']]))
}else{
  print(paste0("make pooled dada object"))
  dada <- dada(derep, err=err, multithread=snakemake@threads,
               BAND_SIZE=snakemake@config[['dada']][['band_size']],
               HOMOPOLYMER_GAP_PENALTY=snakemake@config[['dada']][['homopolymer_gap_penalty']],
               OMEGA_A=as.numeric(snakemake@config[['dada']][['omega_A']]),
               OMEGA_P=as.numeric(snakemake@config[['dada']][['omega_P']]),
               OMEGA_C=as.numeric(snakemake@config[['dada']][['omega_C']]),
               GAPLESS=as.logical(snakemake@config[['dada']][['gapless']]),
               KDIST_CUTOFF=as.numeric(snakemake@config[['dada']][['kdist_cutoff']]),
               MATCH=as.numeric(snakemake@config[['dada']][['match']]),
               MISMATCH=as.numeric(snakemake@config[['dada']][['mismatch']]),
               GAP_PENALTY=as.numeric(snakemake@config[['dada']][['gap_penalty']]),
               selfConsist=as.logical(snakemake@config[['dada']][['selfConsist']]),
               pool=TRUE,
               priors=priors,
               errorEstimationFunction=match.fun(snakemake@config[['dada']][['errorEstimationFunction']]))
  
}
saveRDS(dada,mergefile)
print("done")
