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
library(dada2)
library(Biostrings)


print("merging runs")
# # Merge multiple runs (if necessary)
seqTabList <- unlist(snakemake@input)
print(length(seqTabList))
if(length(seqTabList)>1){
  seqTabs <- list()
  for(i in 1:length(seqTabList)){
     tmpTab <- readRDS(seqTabList[i])
     tabname <- gsub("sequenceTables/seqTab.","",gsub(".RDS","",seqTabList[i]))
     if(all(dim(tmpTab)>0)){
       seqTabs[[tabname]] <- tmpTab
     }else{
       print(paste0("No reads were recovered from run ",tabname,"!")) 
     }
  }
  if(length(seqTabs)>1){
    seqtab <- mergeSequenceTables(tables=seqTabs,repeats="sum")
  }else{
    seqtab <- seqTabs[[1]]
  }
}else{
  seqtab <- readRDS(seqTabList[1])
}
if(snakemake@config[['chimeras']][['remove']]){
  print("Removing chimeras")
  saveRDS(seqtab,snakemake@output[[4]])
  seqs <- DNAStringSet(colnames(seqtab))
  names(seqs) <- sprintf("OTU_%06d",1:length(seqs))
  writeXStringSet(seqs,snakemake@output[[5]])
  write.table(t(seqtab),snakemake@output[[6]],col.names=NA,sep="\t",quote=F)
  seqtab <- removeBimeraDenovo(seqtab, method=snakemake@config[['chimeras']][['method']])
}
print("Saving merged OTU table")
saveRDS(seqtab,snakemake@output[[1]])
seqs <- DNAStringSet(colnames(seqtab))
names(seqs) <- sprintf("OTU_%06d",1:length(seqs))
writeXStringSet(seqs,snakemake@output[[2]])
write.table(t(seqtab),snakemake@output[[3]],col.names=NA,sep="\t",quote=F)


