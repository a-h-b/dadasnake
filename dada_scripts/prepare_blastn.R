log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

.libPaths(paste0(snakemake@config[['dada_lib']],"/R/library"))

library(Biostrings)


print("reading OTU table")
seqTab <- readRDS(snakemake@input[[1]])
if(snakemake@config[['taxonomy']][['decipher']][['do']] &  snakemake@config[['taxonomy']][['mothur']][['do']]){
  missingTax <- DNAStringSet(seqTab$Row.names[seqTab$domain.mothur =="" & seqTab$domain.decipher == ""])
  names(missingTax) <- seqTab$OTU[seqTab$domain.mothur =="" & seqTab$domain.decipher == ""]
}else if(!snakemake@config[['taxonomy']][['decipher']][['do']] &  snakemake@config[['taxonomy']][['mothur']][['do']]){
  missingTax <- DNAStringSet(seqTab$Row.names[seqTab$domain.mothur ==""])
  names(missingTax) <- seqTab$OTU[seqTab$domain.mothur ==""]
}else if(snakemake@config[['taxonomy']][['decipher']][['do']] &  !snakemake@config[['taxonomy']][['mothur']][['do']]){
  missingTax <- DNAStringSet(seqTab$Row.names[seqTab$domain.decipher == ""])
  names(missingTax) <- seqTab$OTU[seqTab$domain.decipher == ""]
}else{
  missingTax <- DNAStringSet(seqTab$Row.names)
  names(missingTax) <- seqTab$OTU
}
print("Saving OTUs for blastn")
writeXStringSet(missingTax,snakemake@output[[1]])


