log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(Biostrings)


print("reading OTU table")
seqTab <- readRDS(snakemake@input[[1]])
print("filtering")
filtCols <- which(grepl("^Domain",colnames(seqTab))|grepl("^Level_1",colnames(seqTab)))
miss <- which(apply(seqTab[,filtCols],1,function(x) all(is.na(x)|x=="")))
missingTax <- DNAStringSet(seqTab$Row.names[miss])
names(missingTax) <- seqTab$OTU[miss]
print("Saving OTUs for blastn")
writeXStringSet(missingTax,snakemake@output[[1]])


