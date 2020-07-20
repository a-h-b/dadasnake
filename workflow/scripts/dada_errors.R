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

filts <- unlist(snakemake@input)

sizes <- sapply(filts,function(x) as.numeric(unlist(strsplit(system2("zcat",args=c(x,"| wc -l"),stdout=T),split=" "))[1]))
filts <- filts[sizes>0]

errfile <- snakemake@output[[1]]
errpath <- gsub("/.+","",errfile) 

set.seed(snakemake@config[['error_seed']])

print(paste0("learning error models ",snakemake@wildcards[['direction']]))
if(length(filts)>0){
    errs <- learnErrors(filts, multithread=snakemake@threads)
    saveRDS(errs,errfile)
pdf(snakemake@output[[2]],width=8,height=11,pointsize=7)
    plotErrors(errs, nominalQ=TRUE)
dev.off()
}else{
   stop("No reads in any of the input files for learning errors.")
}

print("done")
