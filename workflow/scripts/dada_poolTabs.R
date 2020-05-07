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
library(Biostrings)

dadaFile <- snakemake@input[[1]]
dadalist <- readRDS(dadaFile)
seqtab1 <- makeSequenceTable(dadalist)
if(nrow(seqtab1)>1){
  seqtab2 <- aggregate(seqtab1,list(gsub(".+/","",rownames(seqtab1))),sum)
  rm(seqtab1)
  seqtab <- as.matrix(seqtab2[,-1])
  rownames(seqtab) <- seqtab2[,1]
  rm(seqtab2)
}else{
  seqtab <- seqtab1
  rownames(seqtab) <- gsub(".+/","",rownames(seqtab1))
  rm(seqtab1)
}

if(snakemake@config[['chimeras']][['remove']]){
  print("Removing chimeras")
  saveRDS(seqtab,snakemake@output[[5]])
  seqs1 <- DNAStringSet(colnames(seqtab))
  seqtab1 <- seqtab
  seqtab <- removeBimeraDenovo(seqtab, 
                               method=snakemake@config[['chimeras']][['method']], 
                               minFoldParentOverAbundance=snakemake@config[['chimeras']][['minFoldParentOverAbundance']], 
                               minParentAbundance=snakemake@config[['chimeras']][['minParentAbundance']], 
                               allowOneOff=as.logical(snakemake@config[['chimeras']][['allowOneOff']]), 
                               minOneOffParentDistance=snakemake@config[['chimeras']][['minOneOffParentDistance']], 
                               maxShift=snakemake@config[['chimeras']][['maxShift']])
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
writeXStringSet(seqs,snakemake@output[[3]])
outtab <- merge(t(seqtab),data.frame("OTU"=names(seqs),"seq"=seqs,stringsAsFactors=F),
                   by.x=0,by.y="seq")
saveRDS(outtab,snakemake@output[[2]])
write.table(outtab,snakemake@output[[4]],row.names=F,sep="\t",quote=F)

