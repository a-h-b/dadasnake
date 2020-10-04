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

merges <- sort(unlist(snakemake@input))
sizes <- sapply(merges,function(x) file.info(x)$size)
merges <- merges[sizes>0]

sampleNames <- gsub("/","",gsub(".+/","",gsub(".RDS","",merges)))

dadafile <- snakemake@output[[1]]
dadatab <- snakemake@output[[2]]
dadatabtsv <- snakemake@output[[3]]

# Merge
print("merge dada objects")
mergers <- vector("list", length(sampleNames))
names(mergers) <- sampleNames
for(sam in merges) {
  print(paste0("Processing: ", sam))
  csam <- gsub("/","",gsub(".+/","",gsub(".RDS","",sam)))
  mergers[[csam]] <- readRDS(sam)
}
saveRDS(mergers,dadafile)
print("making sequence tab")
# Construct sequence table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, dadatab)
write.table(t(seqtab),dadatabtsv,col.names=NA,sep="\t",quote=F)

print("done")
