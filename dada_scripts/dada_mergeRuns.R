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
  saveRDS(seqtab,snakemake@output[[5]])
  seqs1 <- DNAStringSet(colnames(seqtab))
  seqtab1 <- seqtab
#  names(seqs1) <- sprintf("OTU_%06d",1:length(seqs1))
#  write.table(t(seqtab),snakemake@output[[6]],col.names=NA,sep="\t",quote=F)
  seqtab <- removeBimeraDenovo(seqtab, method=snakemake@config[['chimeras']][['method']])
  seqs <- DNAStringSet(colnames(seqtab))
  names(seqs) <- sprintf("OTU_%06d",1:length(seqs))
  seqs1set <- append(seqs,setdiff(seqs1,seqs))
  names(seqs1set)[names(seqs1set)==""] <- sprintf("chimeric_OTU_%06d",(length(seqs)+1):length(seqs1set))
  writeXStringSet(seqs1set,snakemake@output[[6]])
  outtab1 <- merge(t(seqtab1),data.frame("OTU"=names(seqs1set),"seq"=seqs1set,stringsAsFactors=F),
                   by.x=0,by.y="seq")
  write.table(outtab1,snakemake@output[[7]],row.names=F,sep="\t",quote=F)
}else{
  seqs <- DNAStringSet(colnames(seqtab))
  names(seqs) <- sprintf("OTU_%06d",1:length(seqs))
}
print("Saving merged OTU table")
saveRDS(seqtab,snakemake@output[[1]])
#seqs <- DNAStringSet(colnames(seqtab))
#names(seqs) <- sprintf("OTU_%06d",1:length(seqs))
writeXStringSet(seqs,snakemake@output[[3]])
outtab <- merge(t(seqtab),data.frame("OTU"=names(seqs),"seq"=seqs,stringsAsFactors=F),
                   by.x=0,by.y="seq")
saveRDS(outtab,snakemake@output[[2]])
write.table(outtab,snakemake@output[[4]],row.names=F,sep="\t",quote=F)


