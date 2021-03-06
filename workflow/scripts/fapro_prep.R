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

# File parsing
seqTab <- readRDS(snakemake@input[[1]])
sInfo <- read.delim(snakemake@input[[2]],stringsAsFactors=F,row.names=1)

seqMat <- seqTab[,colnames(seqTab) %in% c("Row.names",rownames(sInfo))]
seqMat$taxonomy <- seqTab[,paste0("taxonomy.",snakemake@config[["postprocessing"]][["faprotax"]][["classifier"]])]
seqMat$OTU <- seqTab$OTU
write.table(seqMat,snakemake@output[[1]],row.names=F,sep="\t",quote=F)


