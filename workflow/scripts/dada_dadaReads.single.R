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
if(!require(dada2)){
  BiocManager::install("GenomeInfoDbData",update=F,ask=F)
  require(dada2)
}
source(snakemake@params[['errorFunctions']])

# File parsing
errfile <- snakemake@input[[1]]

filt <- snakemake@input[[2]]

sampleName <- gsub(".+/","",gsub(".fastq.gz","",filt))

names(filt) <- sampleName

mergefile <- snakemake@output[[1]]

if(snakemake@config[['dada']][['priors']]!=""){
  priorFasta <- readDNAStringSet(snakemake@config[['dada']][['priors']])
  priors <- as.character(priorFasta)
}else{
  priors <- ""
}

if(as.logical(snakemake@config[['dada']][['no_error_assumptions']])){
  err <- NULL
}else{
  err <- readRDS(errfile)
}


# Sample inference and merger of paired-end reads
if(as.numeric(unlist(strsplit(system2("zcat",args=c(filt," 2>/dev/null | head | wc -l"),stdout=T),split=" "))[1])>4){
print(paste0("make dada object, ",sampleName))
derep <- derepFastq(filt)
if("list" %in% class(derep)){
  derep <- lapply(derep, function(x){
    if(any(x$quals[!is.na(x$quals)]<0)){ x$quals <- x$quals+31 }
    return(x)
  })
} else {
  if(any(derep$quals[!is.na(derep$quals)]<0)) derep$quals <- derep$quals+31
}
dada <- dada(derep, err=err, multithread=snakemake@threads,
             BAND_SIZE=as.numeric(snakemake@config[['dada']][['band_size']]),
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
             pool=FALSE,
             priors=priors,
             errorEstimationFunction=match.fun(snakemake@config[['dada']][['errorEstimationFunction']]),
             USE_QUALS=as.logical(snakemake@config[['dada']][['use_quals']]))
saveRDS(dada,mergefile)
}else{
  system(paste("touch", mergefile))
}
print("done")
