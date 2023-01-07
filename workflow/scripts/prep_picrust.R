log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))


print("reading ASV table")
seqTab <- readRDS(snakemake@input[[1]])
if(any(colnames(seqTab)=="ASV")){
 colnames(seqTab)[which(colnames(seqTab)=="ASV")] <- "#OTU_ID"
}else{
 seqTab$OTU <- gsub("OTU_","ASV_",seqTab$OTU)
 colnames(seqTab)[which(colnames(seqTab)=="OTU")] <- "#OTU_ID"
}

print("write ASV table for picrust")
write.table(seqTab[,c("#OTU_ID",colnames(seqTab)[which(sapply(seqTab,class) %in% c("integer","numeric"))])],
  snakemake@output[[1]],sep="\t",quote=F,row.names=F)

