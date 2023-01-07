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

print("filtering")
filtCols <- which(grepl("^Domain",colnames(seqTab))|grepl("^Level_1",colnames(seqTab)))
miss <- which(apply(seqTab[,filtCols],1,function(x) all(is.na(x)|x=="")))
missingTax <- DNAStringSet(seqTab$Row.names[miss])
names(missingTax) <- seqTab[miss,which(colnames(seqTab)==units)]
print("Saving sequences for blastn")
writeXStringSet(missingTax,snakemake@output[[1]])


