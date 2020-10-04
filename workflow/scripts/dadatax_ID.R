log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(BiocParallel)
if (snakemake@threads > 1) {
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}else{
    parallel <- FALSE
    register(SerialParam())
}
library(DECIPHER)
library(Biostrings)
if(!require(dada2)){
  BiocManager::install("GenomeInfoDbData",update=F,ask=F)
  require(dada2)
}

print("finding training set:")
theoPath <- paste0(snakemake@config[['taxonomy']][['dada']][['db_path']],"/",snakemake@config[['taxonomy']][['dada']][['refFasta']])
if(file.exists(theoPath)){
  loadPath <- theoPath
}else{
  altPath <- paste0(snakemake@config[['taxonomy']][['decipher']][['db_path']],"/",
                    list.files(path=snakemake@config[['taxonomy']][['dada']][['db_path']],
                               pattern=snakemake@config[['taxonomy']][['dada']][['refFasta']]))
  if(length(altPath)==1) loadPath <- altPath else stop(paste0(theoPath," not found."))
}
print(loadPath)
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus")

print(paste0("loading sequence ",snakemake@input))
seqs <- readDNAStringSet(snakemake@input[[1]])

print("running classification")
set.seed(snakemake@config[['taxonomy']][['dada']][['seed']])
ids <- assignTaxonomy(as.character(seqs), loadPath, tryRC=snakemake@config[['taxonomy']][['dada']][['tryRC']],
               minBoot=snakemake@config[['taxonomy']][['dada']][['minBoot']], 
               multithread=parallel, taxLevels=ranks, verbose=T)
if(snakemake@config[['taxonomy']][['dada']][['look_for_species']]){
  ids <- addSpecies(ids,snakemake@config[['taxonomy']][['dada']][['spec_db']])
}
rownames(ids) <- names(seqs)
ids <- ids[order(rownames(ids)),]

write.table(data.frame("taxonomy"=apply(ids,1,function(x) paste(x[!is.na(x)],sep=";",collapse=";")),
                       ids,"seq"=as.character(seqs)[order(names(seqs))],stringsAsFactors=F),
            snakemake@output[[1]],sep="\t",col.names=NA,quote=F)
saveRDS(data.frame("taxonomy"=apply(ids,1,function(x) paste(x[!is.na(x)],sep=";",collapse=";")),
                       ids,"seq"=as.character(seqs)[order(names(seqs))],stringsAsFactors=F),
        snakemake@output[[2]])

