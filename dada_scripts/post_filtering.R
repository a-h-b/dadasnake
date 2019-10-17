log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

.libPaths(paste0(snakemake@config[['dada_lib']],"/R/library"))

library(BiocParallel)
#parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}else{
    parallel <- FALSE
    register(SerialParam())
}
#register(SerialParam())
library(Biostrings)

# File parsing
seqs <- readDNAStringSet((snakemake@input[[1]])
seqTab <- readRDS(snakemake@input[[2]])

if(snakemake@config[['final_table_filtering']][['keep_target_taxa']]!="."){
 if(any(grepl("taxonomy",colnames(seqTab)))){
  taxcols <- grep("taxonomy",colnames(seqTab))
  if(length(taxcols)>1){
   taxS <- apply(seqTab[,taxcols],1,function(x)paste(x,collapse=";",sep=";"))
  }else{
   taxS <- seqTab[,taxcols]
  }
  keep <- grep(snakemake@config[['final_table_filtering']][['keep_target_taxa']],taxS)
  filtTab <- seqTab[keep,]
 }else{
  print(paste0("No taxonomy in sequence table ",snakemake@input[[2]],"."))
  filtTab <- seqTab
 }
}else{
 filtTab <- seqTab
}
if(nrow(filtTab)==0){
 print("No OTUs left after taxonomy filter.")
 system(paste0("touch ",snakemake@output[[1]]))
 system(paste0("touch ",snakemake@output[[2]]))
 system(paste0("touch ",snakemake@output[[3]]))
}else{
 slen <- nchar(filtTab$Row.names)
 keep <- which(slen >= snakemake@config[['final_table_filtering']][['target_min_length']] 
               & slen <= snakemake@config[['final_table_filtering']][['target_max_length']])
 filtTab <- filtTab[keep,]
 if(nrow(filtTab)==0){
  print("No OTUs left after length filter.")
  system(paste0("touch ",snakemake@output[[1]]))
  system(paste0("touch ",snakemake@output[[2]]))
  system(paste0("touch ",snakemake@output[[3]]))
 }else{
  writeRDS(filtTab,snakemake@output[[1]])
  write.table(filtTab,snakemake@output[[2]],row.names=F,sep="\t",quote=F)
  filtSeq <- seqs[[filtTab$OTU]]
  writeXStringSet(filtSeq,snakemake@output[[3]])
 }
}

