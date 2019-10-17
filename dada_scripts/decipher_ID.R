log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

.libPaths(paste0(snakemake@config[['dada_lib']],"/R/library"))

library(BiocParallel)
#parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}else{
    parallel <- FALSE
    register(SerialParam())
}
#register(SerialParam())
library(DECIPHER)
library(Biostrings)

print("loading training set:")
theoPath <- paste0(snakemake@config[['taxonomy']][['decipher']][['db_path']],"/",snakemake@config[['taxonomy']][['decipher']][['tax_db']])
if(file.exists(theoPath)){
  loadPath <- theoPath
}else{
  altPath <- paste0(snakemake@config[['taxonomy']][['decipher']][['db_path']],"/",
                    list.files(path=snakemake@config[['taxonomy']][['decipher']][['db_path']],
                               pattern=snakemake@config[['taxonomy']][['decipher']][['tax_db']]))
  if(length(altPath)==1) loadPath <- altPath else stop(paste0(theoPath," not found."))
}
print(loadPath)
load(loadPath)
ranks <- c("Domain","Phylum","Class","Order","Family","Genus")

print(paste0("loading sequence ",snakemake@input))
seqs <- readDNAStringSet(snakemake@input[[1]])

print("running classification")
set.seed(snakemake@config[['taxonomy']][['decipher']][['seed']])
ids <- IdTaxa(seqs, trainingSet, strand=snakemake@config[['taxonomy']][['decipher']][['strand']],
               threshold=snakemake@config[['taxonomy']][['decipher']][['threshold']],
               bootstraps=snakemake@config[['taxonomy']][['decipher']][['bootstraps']], 
               processors=snakemake@threads, verbose=T)
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- names(seqs)
if(snakemake@config[['taxonomy']][['decipher']][['look_for_species']]){
  taxid <- addSpecies(taxid,snakemake@config[['taxonomy']][['decipher']][['spec_db']])
}
taxid <- taxid[order(rownames(taxid)),]

write.table(data.frame("taxonomy"=apply(taxid,1,function(x) paste(x[!is.na(x)],sep=";",collapse=";")),
                       taxid,"seq"=as.character(seqs)[order(names(seqs))],stringsAsFactors=F),
            snakemake@output[[1]],sep="\t",col.names=NA,quote=F)
saveRDS(data.frame("taxonomy"=apply(taxid,1,function(x) paste(x[!is.na(x)],sep=";",collapse=";")),
                       taxid,"seq"=as.character(seqs)[order(names(seqs))],stringsAsFactors=F),
        snakemake@output[[2]])

