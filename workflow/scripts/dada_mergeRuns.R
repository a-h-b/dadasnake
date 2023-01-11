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
library(Biostrings)


print("merging runs")
# # Merge multiple runs (if necessary)
seqTabList <- unlist(snakemake@input)
print(length(seqTabList))
if(length(seqTabList)>1){
  seqTabs <- list()
  for(i in 1:length(seqTabList)){
     if(file.info(seqTabList[i])$size > 0){
       tmpTab <- readRDS(seqTabList[i])
       tabname <- gsub("sequenceTables/seqTab.","",gsub(".RDS","",seqTabList[i]))
       print(tabname)
       if(all(dim(tmpTab)>0)){
         seqTabs[[tabname]] <- tmpTab
       }else{
         print(paste0("No reads were recovered from run ",tabname,"!")) 
       }
     }
  }
  if(length(seqTabs)>1){
    print("merging tabs")
    seqtab <- mergeSequenceTables(tables=seqTabs,repeats="sum")
  }else{
    seqtab <- seqTabs[[1]]
  }
  rm(seqTabs)
}else{
  seqtab <- readRDS(seqTabList[1])
}
num_digits <- as.character(ceiling(log10(ncol(seqtab))))
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
  names(seqs) <- sprintf(paste0("ASV_%0",num_digits,"d"),1:length(seqs))
  seqs1set <- append(seqs,setdiff(seqs1,seqs))
  names(seqs1set)[names(seqs1set)==""] <- sprintf(paste0("chimeric_ASV_%0",num_digits,"d"),(length(seqs)+1):length(seqs1set))
  writeXStringSet(seqs1set,snakemake@output[[6]])
  outtab1 <- merge(t(seqtab1),data.frame("ASV"=names(seqs1set),"seq"=seqs1set,stringsAsFactors=F),
                   by.x=0,by.y="seq")
  write.table(outtab1,snakemake@output[[7]],row.names=F,sep="\t",quote=F)
  rm(outtab1)
  rm(seqs1set)
  rm(seqtab1)
}else{
  seqs <- DNAStringSet(colnames(seqtab))
  names(seqs) <- sprintf(paste0("ASV_%0",num_digits,"d"),1:length(seqs))
}
print("Saving sequences")
saveRDS(seqtab,snakemake@output[[1]])
writeXStringSet(seqs,snakemake@output[[3]])
print("Saving merged ASV table")
sams <- rownames(seqtab)
seqtab <- data.frame(t(seqtab),stringsAsFactors=F)
colnames(seqtab) <- sams
seqtab$Row.names <- rownames(seqtab)
seqtab$ASV <- names(seqs)
saveRDS(seqtab[,c(ncol(seqtab)-1,1:(ncol(seqtab)-2),ncol(seqtab))],snakemake@output[[2]])
write.table(seqtab[,c(ncol(seqtab)-1,1:(ncol(seqtab)-2),ncol(seqtab))],snakemake@output[[4]],row.names=F,sep="\t",quote=F)


