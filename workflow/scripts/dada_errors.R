log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(BiocParallel)
if (snakemake@threads > 1) {
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

filts <- unlist(snakemake@input)

sizes <- sapply(filts,function(x) as.numeric(unlist(strsplit(system2("zcat",args=c(x," 2>/dev/null | head |  wc -l"),stdout=T),split=" "))[1]))
filts <- filts[sizes>0]

errfile <- snakemake@output[[1]]
errpath <- gsub("/.+","",errfile) 

set.seed(snakemake@config[['error_seed']])

print(paste0("learning error models ",snakemake@wildcards[['direction']]))
if(length(filts)>0){
    #dereps <- derepFastq(filts)
    NBASES <- 0
    NREADS <- 0
    drps <- vector("list", length(filts))
    filts <- sample(filts)
    for (i in seq_along(filts)) {
      print(i)
      print("derepFastq:")
      drps[[i]] <- derepFastq(filts[[i]])
      if(any( drps[[i]]$quals[!is.na(drps[[i]]$quals)]<0)) drps[[i]]$quals <- drps[[i]]$quals+31
      NREADS <- NREADS + sum(drps[[i]]$uniques)
      NBASES <- NBASES + sum(drps[[i]]$uniques * nchar(names(drps[[i]]$uniques)))
      if (NBASES > as.numeric(snakemake@config[['dada']][['error_nbases']])) {
        break
      }
    }
    drps <- drps[1:i]
    print("learnErrors:")
    errs <- learnErrors(drps, multithread=snakemake@threads, 
                         nbases=as.numeric(snakemake@config[['dada']][['error_nbases']]),
                         MAX_CONSIST=as.numeric(snakemake@config[['dada']][['error_max_consist']]),
                         OMEGA_C=as.numeric(snakemake@config[['dada']][['error_omega_C']]),
                         errorEstimationFunction=match.fun(snakemake@config[['dada']][['errorEstimationFunction']]))
    saveRDS(errs,errfile)
    pdf(snakemake@output[[2]],width=8,height=11,pointsize=7)
       print(plotErrors(errs, nominalQ=TRUE))
    dev.off()
}else{
   stop("No reads in any of the input files for learning errors.")
}

print("done")
