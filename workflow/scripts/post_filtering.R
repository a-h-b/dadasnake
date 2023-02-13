log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(BiocParallel)
#parallel <- FALSE
if (snakemake@threads > 1) {
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}else{
    parallel <- FALSE
    register(SerialParam())
}
#register(SerialParam())
library(Biostrings)


if(snakemake@params[["what"]]=="ASV"){
  units <- "ASV"
}else{
  units <- "OTU"
}

# File parsing
seqs <- readDNAStringSet(snakemake@input[[1]])
if(units=="ASV" & any(grepl("^OTU",names(seqs)))){
  names(seqs) <- gsub("OTU_","ASV_",names(seqs))
}
seqTab <- readRDS(snakemake@input[[2]])
if(units=="ASV" & any(colnames(seqTab)=="OTU")){
  seqTab$OTU <- gsub("OTU_","ASV_",seqTab$OTU)
  colnames(seqTab)[which(colnames(seqTab)=="OTU")] <- "ASV"
}
 
if(!snakemake@params[['target_taxa']] %in% c(".","")){
 if(any(grepl("^taxonomy",colnames(seqTab)))){
  taxcols <- grep("^taxonomy",colnames(seqTab))
  if(length(taxcols)>1){
   taxS <- apply(seqTab[,taxcols],1,function(x)paste(x,collapse=";",sep=";"))
  }else{
   taxS <- seqTab[,taxcols]
  }
  keep <- grep(snakemake@params[['target_taxa']],taxS)
  filtTab <- seqTab[keep,]
 }else{
  print(paste0("No taxonomy in sequence table ",snakemake@input[[2]],"."))
  filtTab <- seqTab
 }
}else{
 filtTab <- seqTab
}
if(nrow(filtTab)==0){
 print("No entries left after taxonomy filter.")
 system(paste0("touch ",snakemake@output[[1]]))
 system(paste0("touch ",snakemake@output[[2]]))
 system(paste0("touch ",snakemake@output[[3]]))
}else{
 slen <- nchar(filtTab$Row.names)
 keep <- which(slen >= snakemake@params[['min']] 
               & slen <= snakemake@params[['max']])
 filtTab <- filtTab[keep,]
 if(nrow(filtTab)==0){
  print("No entries left after length filter.")
  system(paste0("touch ",snakemake@output[[1]]))
  system(paste0("touch ",snakemake@output[[2]]))
  system(paste0("touch ",snakemake@output[[3]]))
 }else{
  saveRDS(filtTab,snakemake@output[[1]])
  write.table(filtTab,snakemake@output[[2]],row.names=F,sep="\t",quote=F)
  filtSeq <- seqs[filtTab[,which(colnames(filtTab)==units)]]
  writeXStringSet(filtSeq,snakemake@output[[3]])
 }
}

