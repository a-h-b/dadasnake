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
print(paste(length(seqTabList), "runs to read"))

# if there are more than 1 run, a list for pre (and post) chimera-removal ASV is made) 
if(length(seqTabList)>1){
  seqTabs <- list()
  if(snakemake@config[['chimeras']][['remove_by_run']]){
    seqTabs_prechim <- list()
  }
  # read table for each run and remove chimera, if applicable, by run
  for(i in 1:length(seqTabList)){
     if(file.info(seqTabList[i])$size > 0){
       tmpTab <- readRDS(seqTabList[i])
       tabname <- gsub("sequenceTables/seqTab.","",gsub(".RDS","",seqTabList[i]))
       print(tabname)
       if(all(dim(tmpTab)>0)){
         seqTabs[[tabname]] <- tmpTab
         if(snakemake@config[['chimeras']][['remove_by_run']]){
           seqTabs_prechim[[tabname]] <- tmpTab
           print(paste("Removing chimeras for",tabname))
           seqTabs[[tabname]] <- removeBimeraDenovo(tmpTab,
                                        method=snakemake@config[['chimeras']][['method']],
                                        minFoldParentOverAbundance=snakemake@config[['chimeras']][['minFoldParentOverAbundance']],
                                        minParentAbundance=snakemake@config[['chimeras']][['minParentAbundance']],
                                        allowOneOff=as.logical(snakemake@config[['chimeras']][['allowOneOff']]),
                                        minOneOffParentDistance=snakemake@config[['chimeras']][['minOneOffParentDistance']],
                                        maxShift=snakemake@config[['chimeras']][['maxShift']])
         }
       }else{
         print(paste0("No reads were recovered from run ",tabname,"!")) 
       }
     }
  }
  if(length(seqTabs)>1){
    # merge final tables
    print("merging tabs")
    seqtab <- mergeSequenceTables(tables=seqTabs,repeats="sum")
  }else{
    seqtab <- seqTabs[[1]]
  }
  rm(seqTabs)
  if(snakemake@config[['chimeras']][['remove_by_run']] & length(seqTabs_prechim)>1){
    # merge pre-chimera tabs, if applicable
    print("merging prechimera tabs")
    seqtab_prechim <- mergeSequenceTables(tables=seqTabs_prechim,repeats="sum")
    saveRDS(seqtab_prechim,snakemake@output[[5]])
  }else{
    if(snakemake@config[['chimeras']][['remove_by_run']]){
      seqtab_prechim <- seqTabs_prechim[[1]]
      saveRDS(seqtab_prechim,snakemake@output[[5]])
    }
  }
}else{
  # if there's just one input, this is used as table
  if(file.info(seqTabList[1])$size > 0){
    seqtab <- readRDS(seqTabList[1])
    if(snakemake@config[['chimeras']][['remove']]){
      saveRDS(seqtab,snakemake@output[[5]])
      seqtab_prechim <- seqtab
    }
  }
}
# at this point, there's a seqtab - the chimeras might still need to be removed
num_digits <- as.character(ceiling(log10(ncol(seqtab))))
if(snakemake@config[['chimeras']][['remove']] & (length(seqTabList)==1 | !snakemake@config[['chimeras']][['remove_by_run']])){
  print("Removing chimeras")
  seqtab <- removeBimeraDenovo(seqtab, 
                               method=snakemake@config[['chimeras']][['method']],
                               minFoldParentOverAbundance=snakemake@config[['chimeras']][['minFoldParentOverAbundance']],
                               minParentAbundance=snakemake@config[['chimeras']][['minParentAbundance']],
                               allowOneOff=as.logical(snakemake@config[['chimeras']][['allowOneOff']]),
                               minOneOffParentDistance=snakemake@config[['chimeras']][['minOneOffParentDistance']],
                               maxShift=snakemake@config[['chimeras']][['maxShift']])
}
seqs <- DNAStringSet(colnames(seqtab))
names(seqs) <- sprintf(paste0("OTU_%0",num_digits,"d"),1:length(seqs))
if(snakemake@config[['chimeras']][['remove']]){
  seqs_prechim <- DNAStringSet(colnames(seqtab_prechim))
  seqs_prechim <- append(seqs,setdiff(seqs_prechim,seqs))
  names(seqs_prechim)[names(seqs_prechim)==""] <- sprintf(paste0("chimeric_OTU_%0",num_digits,"d"),(length(seqs)+1):length(seqs_prechim))
  writeXStringSet(seqs_prechim,snakemake@output[[6]])
  outtab_prechim <- merge(t(seqtab_prechim),
                          data.frame("OTU"=names(seqs_prechim),"seq"=seqs_prechim,stringsAsFactors=F),
                          by.x=0,by.y="seq")
  write.table(outtab_prechim,snakemake@output[[7]],row.names=F,sep="\t",quote=F)
  rm(outtab_prechim)
  rm(seqs_prechim)
  rm(seqtab_prechim)
}
print("Saving sequences")
saveRDS(seqtab,snakemake@output[[1]])
writeXStringSet(seqs,snakemake@output[[3]])
print("Saving merged OTU table")
seqtab <- data.frame(t(seqtab),stringsAsFactors=F)
seqtab$Row.names <- rownames(seqtab)
seqtab$OTU <- names(seqs)
saveRDS(seqtab[,c(ncol(seqtab)-1,1:(ncol(seqtab)-2),ncol(seqtab))],snakemake@output[[2]])
write.table(seqtab[,c(ncol(seqtab)-1,1:(ncol(seqtab)-2),ncol(seqtab))],snakemake@output[[4]],row.names=F,sep="\t",quote=F)


