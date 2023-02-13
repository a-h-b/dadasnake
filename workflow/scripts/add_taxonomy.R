log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(Biostrings)


if(snakemake@params[["what"]]=="ASV"){
  units <- "ASV"
}else{
  units <- "OTU"
}

print("reading OTU table")
seqTab <- readRDS(snakemake@input[[1]])
if(units=="ASV" & any(colnames(seqTab)=="OTU")){
  seqTab$OTU <- gsub("OTU_","ASV_",seqTab$OTU)
  colnames(seqTab)[which(colnames(seqTab)=="OTU")] <- "ASV"
}

if(snakemake@params[['ITSx']]=="add"){
  print("reading ITSx result")
  iadd <- 1
  if(as.numeric(file.info(snakemake@input[[2]])$size)>0){
    seqs <- readDNAStringSet(snakemake@input[[2]])
    if(units=="ASV" & any(grepl("^OTU",names(seqs)))){
      names(seqs) <- gsub("OTU_","ASV_",names(seqs))
    }
    seqTab$ITSx <- seqTab[,which(colnames(seqTab)==units)] %in% names(seqs)
    rm(seqs)
  }
}else{
  iadd <- 0
}

if(length(snakemake@input)>=2+iadd){
  print("reading taxonomy results")
  for(i in (2+iadd):length(snakemake@input)){
    cTax <- readRDS(snakemake@input[[i]])
    if(units=="ASV"){
      if(any(colnames(cTax)=="OTU")){
        cTax$OTU <- gsub("OTU_","ASV_",cTax$OTU)
        colnames(cTax)[which(colnames(cTax)=="OTU")] <- "ASV"
      }else{
        if(any(grepl("^OTU",cTax$ASV))) cTax$ASV <- gsub("OTU_","ASV_",cTax$ASV)
      }
    }
    seqTab <- merge(seqTab,cTax,by=units,all.x=T) 
  }
}
seqTab[is.na(seqTab)] <- ""
print("Saving OTU table with taxonomy")
saveRDS(seqTab,snakemake@output[[2]])
write.table(seqTab,snakemake@output[[1]],row.names=F,sep="\t",quote=F)


