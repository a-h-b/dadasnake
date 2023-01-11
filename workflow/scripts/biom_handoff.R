log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(BiocParallel)
if (snakemake@threads > 1) {
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}else{
    parallel <- FALSE
    register(SerialParam())
}
library(biomformat)


# File parsing
seqTab <- readRDS(snakemake@input[[1]])
if(any(colnames(seqTab)=="OTU")){
  seqTab$OTU <- gsub("OTU_","ASV_",seqTab$OTU)
  colnames(seqTab)[which(colnames(seqTab)=="OTU")] <- "ASV"
}
sInfo <- read.delim(snakemake@input[[2]],stringsAsFactors=F,row.names=1)
if(!any(grepl("-",colnames(seqTab)))) rownames(sInfo) <- gsub("-",".",rownames(sInfo))
if(length(which(colnames(seqTab) %in% ifelse(grepl("^[[:digit:]]",rownames(sInfo)),
                                                             paste0("X",rownames(sInfo)),
                                                              rownames(sInfo)))) > 0){
  if(snakemake@params[["currentStep"]]=="dada"){
   seqMat <- seqTab[,-c(1,ncol(seqTab))]
   rownames(seqMat) <- seqTab[,1]
   seqMeta <- data.frame("OTU_ID"=seqTab$ASV,stringsAsFactors = F,row.names = seqTab[,1])
  }else if(snakemake@params[["currentStep"]]=="taxonomy"){
   seqMat <- seqTab[,colnames(seqTab) %in% ifelse(grepl("^[[:digit:]]",rownames(sInfo)),
                                                             paste0("X",rownames(sInfo)),
                                                              rownames(sInfo))]
   if(class(seqMat) == "integer"){
     seqMat <- matrix(seqMat,ncol=1)
     colnames(seqMat) <- intersect(colnames(seqTab), ifelse(grepl("^[[:digit:]]",rownames(sInfo)),
                                                             paste0("X",rownames(sInfo)),
                                                              rownames(sInfo)))
   }
   rownames(seqMat) <- seqTab$Row.names
   seqMeta <- seqTab[,!colnames(seqTab) %in% c(ifelse(grepl("^[[:digit:]]",rownames(sInfo)),
                                                             paste0("X",rownames(sInfo)),
                                                              rownames(sInfo)),"Row.names")]
  }
  sInfo <- sInfo[ifelse(grepl("^[[:digit:]]",rownames(sInfo)),
                         paste0("X",rownames(sInfo)),
                         rownames(sInfo)) %in% colnames(seqMat),]
  seqBiom <- make_biom(seqMat,observation_metadata = seqMeta,sample_metadata=sInfo)
  write_biom(seqBiom,snakemake@output[[1]])
}else{
  print("Nothing to write.")
  write.table("",snakemake@output[[1]],quote=F,col.names=F,row.names=F)
}

