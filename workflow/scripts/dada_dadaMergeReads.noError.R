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
if(!require(dada2)){
  BiocManager::install("GenomeInfoDbData",update=F,ask=F)
  require(dada2)
}

# File parsing
filtF <- snakemake@input[[1]]
filtR <- snakemake@input[[2]]

sampleNameF <- gsub("/","",gsub("downsampled/","",gsub("filtered/","run.",gsub(".fwd.fastq.gz","",filtF))))
sampleNameR <- gsub("/","",gsub("downsampled/","",gsub("filtered/","run.",gsub(".rvs.fastq.gz","",filtR))))

if(sampleNameF!=sampleNameR) stop("Forward and reverse files do not match.")

sampleName <- gsub(".+/","",gsub(".fwd.fastq.gz","",filtF))

names(filtF) <- sampleName
names(filtR) <- sampleName

mergefile <- snakemake@output[[1]]

if(as.numeric(unlist(strsplit(system2("zcat",args=c(filtF," 2>/dev/null | head | wc -l"),stdout=T),split=" "))[1])>4){

if(as.logical(snakemake@config[['dada']][['no_error_assumptions']])){
  errF <- NULL
  errR <- NULL
}else{
  errF <- inflateErr(tperr1, 3)
  errR <- inflateErr(tperr1, 3)
}

print("merging")

# Sample inference and merger of paired-end reads
print(paste0("make dada object and merge, ",sampleName))
derepF <- derepFastq(filtF)
if("list" %in% class(derepF)){
  derepF <- lapply(derepF, function(x){
    if(any(x$quals[!is.na(x$quals)]<0)){ x$quals <- x$quals+31 }
    return(x)
  })
} else {
  if(any(derepF$quals[!is.na(derepF$quals)]<0)) derepF$quals <- derepF$quals+31
}

dadaF <- dada(derepF, err=errF, multithread=snakemake@threads,
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
               pool=FALSE,
               priors=priors,
               selfConsist=as.logical(snakemake@config[['dada']][['selfConsist']]),
               errorEstimationFunction=match.fun(snakemake@config[['dada']][['errorEstimationFunction']]),
               USE_QUALS=as.logical(snakemake@config[['dada']][['use_quals']]))
derepR <- derepFastq(filtR)
if("list" %in% class(derepR)){
  derepR <- lapply(derepR, function(x){
    if(any(x$quals[!is.na(x$quals)]<0)){ x$quals <- x$quals+31 }
    return(x)  
  })           
} else {       
  if(any(derepR$quals[!is.na(derepR$quals)]<0)) derepR$quals <- derepR$quals+31
}
dadaR <- dada(derepR, err=errR, multithread=snakemake@threads,
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
               pool=FALSE,
               priors=priors,
               selfConsist=as.logical(snakemake@config[['dada']][['selfConsist']]),
               errorEstimationFunction=match.fun(snakemake@config[['dada']][['errorEstimationFunction']]),
               USE_QUALS=as.logical(snakemake@config[['dada']][['use_quals']]))
merger <- mergePairs(dadaF, derepF, dadaR, derepR,
                     minOverlap=snakemake@config[['pair_merging']][['min_overlap']],
                     maxMismatch=snakemake@config[['pair_merging']][['max_mismatch']],
                     justConcatenate=snakemake@config[['pair_merging']][['just_concatenate']],
                     trimOverhang=snakemake@config[['pair_merging']][['trim_overhang']])
saveRDS(merger,mergefile)
}else{
system(paste("touch",mergefile))
}
print("done")
