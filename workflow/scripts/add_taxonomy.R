log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(Biostrings)


print("reading OTU table")
seqTab <- readRDS(snakemake@input[[1]])
if(snakemake@config[['ITSx']][['do']]){
  print("reading ITSx result")
  iadd <- 1
  seqs <- readDNAStringSet(snakemake@input[[2]])
  seqTab$ITSx <- seqTab$OTU %in% names(seqs)
  rm(seqs)
}else{
  iadd <- 0
}
for(i in (2+iadd):length(snakemake@input)){
  cTax <- readRDS(snakemake@input[[i]])
  seqTab <- merge(seqTab,cTax,by="OTU",all.x=T) 
}
seqTab[is.na(seqTab)] <- ""
print("Saving OTU table with taxonomy")
saveRDS(seqTab,snakemake@output[[2]])
write.table(seqTab,snakemake@output[[1]],row.names=F,sep="\t",quote=F)


