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

if(snakemake@params[["what"]]=="ASV"){
  units <- "ASV"
}else{
  units <- "OTU"
}

print("loading training set:")
theoPath <- snakemake@params[['DB']]
if(file.exists(theoPath)){
  loadPath <- theoPath
}else{
  altPath <- paste0(snakemake@params[['db_path']],"/",
                    list.files(path=snakemake@params[['db_path']],
                               pattern=snakemake@params[['tax_db']]))
  if(length(altPath)==1) loadPath <- altPath else stop(paste0(theoPath," not found."))
}
print(loadPath)
load(loadPath)
ranks <- c("Domain","Phylum","Class","Order","Family","Genus")

print(paste0("loading sequence ",snakemake@input))
seqs <- readDNAStringSet(snakemake@input[[1]])
if(units=="ASV" & any(grepl("^OTU",names(seqs)))) names(seqs) <- gsub("OTU_","ASV_",names(seqs))


print("running classification")
set.seed(snakemake@params[['seed']])
ids <- IdTaxa(seqs, trainingSet, strand=snakemake@params[['strand']],
               threshold=snakemake@params[['threshold']],
               bootstraps=snakemake@params[['bootstraps']], 
               processors=snakemake@threads, verbose=T)
taxid <- t(sapply(ids, function(x) {
        m <- match(tolower(ranks), x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- as.character(seqs)
if(snakemake@params[['do_species']]){
  taxid <- addSpecies(taxid,snakemake@params[['SPEC']])
}
rownames(taxid) <- names(seqs)
taxid <- taxid[order(rownames(taxid)),]

write.table(data.frame("taxonomy"=apply(taxid,1,function(x) paste(x[!is.na(x)],sep=";",collapse=";")),
                       taxid,"seq"=as.character(seqs)[order(names(seqs))],stringsAsFactors=F),
            snakemake@output[[1]],sep="\t",col.names=NA,quote=F)
saveRDS(data.frame("taxonomy"=apply(taxid,1,function(x) paste(x[!is.na(x)],sep=";",collapse=";")),
                       taxid,"seq"=as.character(seqs)[order(names(seqs))],stringsAsFactors=F),
        snakemake@output[[2]])

