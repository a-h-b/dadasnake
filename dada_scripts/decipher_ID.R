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
theoPath <- paste0(snakemake@config[['taxonomy']][['db_path']],"/",snakemake@config[['taxonomy']][['tax_db']])
if(file.exists(theoPath)){
  loadPath <- theoPath
}else{
  altPath <- paste0(snakemake@config[['taxonomy']][['db_path']],"/",
                    list.files(path=snakemake@config[['taxonomy']][['db_path']],
                               pattern=snakemake@config[['taxonomy']][['tax_db']]))
  if(length(altPath)==1) loadPath <- altPath else stop(paste0(theoPath," not found."))
}
print(loadPath)
load(loadPath)
if(snakemake@config[['taxonomy']][['look_for_species']]){
  ranks <- c("domain","phylum","class","order","family","genus","species")
}else{
  ranks <- c("domain","phylum","class","order","family","genus")
}

print(paste0("loading seqTab ",snakemake@input))
seqtab <- readRDS(snakemake@input[[1]])
seqs <- DNAStringSet(colnames(seqtab))
names(seqs) <- sprintf("OTU_%06d",1:length(seqs))

print("running classification")
set.seed(snakemake@config[['taxonomy']][['seed']])
ids <- IdTaxa(seqs, trainingSet, strand=snakemake@config[['taxonomy']][['strand']],
               threshold=snakemake@config[['taxonomy']][['threshold']],
               bootstraps=snakemake@config[['taxonomy']][['bootstraps']], 
               processors=snakemake@threads, verbose=T)
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- names(seqs)
taxid <- taxid[order(rownames(taxid)),]
otu <- t(seqtab)
rownames(otu) <- names(seqs)
otu <- otu[order(rownames(otu)),] 
if(!all(rownames(taxid)==rownames(otu))) stop("mismatch between OTU names in taxonomy and OTU table")

write.table(data.frame(otu,"taxonomy"=apply(taxid,1,function(x) paste(x[!is.na(x)],sep=";",collapse=";")),
                       taxid,"seq"=as.character(seqs)[order(names(seqs))],stringsAsFactors=F),
            snakemake@output[[1]],sep="\t",col.names=NA,quote=F)
saveRDS(data.frame(otu,"taxonomy"=apply(taxid,1,function(x) paste(x[!is.na(x)],sep=";",collapse=";")),
                       taxid,"seq"=as.character(seqs)[order(names(seqs))],stringsAsFactors=F),
        snakemake@output[[2]])

