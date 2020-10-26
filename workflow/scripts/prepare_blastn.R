log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(Biostrings)


print("reading OTU table")
seqTab <- readRDS(snakemake@input[[1]])
if(snakemake@config[['taxonomy']][['decipher']][['do']] &  snakemake@config[['taxonomy']][['mothur']][['do']]){
  missingTax <- DNAStringSet(seqTab$Row.names[(seqTab$Domain.mothur ==""|is.na(seqTab$Domain.mothur) )& is.na(seqTab$Domain.decipher)])
  names(missingTax) <- seqTab$OTU[( seqTab$Domain.mothur ==""|is.na(seqTab$Domain.mothur) ) & is.na(seqTab$Domain.decipher)]
}else if(snakemake@config[['taxonomy']][['decipher']][['do']] &  snakemake@config[['taxonomy']][['dada']][['do']]){
  missingTax <- DNAStringSet(seqTab$Row.names[is.na(seqTab$domain.decipher) & is.na(seqTab$Kingdom.dada)])
  names(missingTax) <- seqTab$OTU[is.na(seqTab$Domain.decipher) & is.na(seqTab$Kingdom.dada)]
}else if(!snakemake@config[['taxonomy']][['decipher']][['do']] &  snakemake@config[['taxonomy']][['mothur']][['do']]){
  missingTax <- DNAStringSet(seqTab$Row.names[ seqTab$Domain.mothur ==""|is.na(seqTab$Domain.mothur) ])
  names(missingTax) <- seqTab$OTU[seqTab$Domain.mothur ==""]
}else if(!snakemake@config[['taxonomy']][['decipher']][['do']] &  snakemake@config[['taxonomy']][['dada']][['do']]){
  missingTax <- DNAStringSet(seqTab$Row.names[is.na(seqTab$Kingdom.dada)])
  names(missingTax) <- seqTab$OTU[is.na(seqTab$Kingdom.dada)]
}else if(snakemake@config[['taxonomy']][['decipher']][['do']] & !snakemake@config[['taxonomy']][['dada']][['do']] & !snakemake@config[['taxonomy']][['mothur']][['do']]){
  missingTax <- DNAStringSet(seqTab$Row.names[is.na(seqTab$Domain.decipher)])
  names(missingTax) <- seqTab$OTU[is.na(seqTab$Domain.decipher)]
}else{
  missingTax <- DNAStringSet(seqTab$Row.names)
  names(missingTax) <- seqTab$OTU
}
print("Saving OTUs for blastn")
writeXStringSet(missingTax,snakemake@output[[1]])


