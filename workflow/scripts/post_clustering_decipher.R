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

print(paste0("loading sequence ",snakemake@input[[1]]))
seqs <- readDNAStringSet(snakemake@input[[1]])

print(paste0("loading ASV table",snakemake@input[[2]]))
seqTab <- readRDS(snakemake@input[[2]])

print("running alignment")
aln <- DECIPHER::AlignSeqs(seqs, processors = snakemake@threads)

print("getting distance")
d <- DECIPHER::DistanceMatrix(aln, processors = snakemake@threads)

cutoff <- as.numeric(snakemake@params[["cutoff"]])
if(cutoff > 1) cutoff <- cutoff/100
if(cutoff > 0.5) cutoff <- 1 - cutoff # no one is going to be really wanting to find clusters with less than 50% similarity, so we can catch confused notation 

print(paste("cluster ASVs at", cutoff, "%"))
clusters <- DECIPHER::TreeLine(
  myDistMatrix=d,
  method = "complete",
  cutoff = cutoff, # use `cutoff = 0.03` for a 97% OTU
  type = "clusters",
  processors = snakemake@threads)
colnames(clusters) <- "decipherCluster"


write.table(cluster,
            snakemake@output[[1]],sep="\t",col.names=NA,quote=F)

print("calculate clustered table")
seqTab <- merge(clusters,seqTab,by.x=0,by.y="OTU")

clusTab <- aggregate(seqTab[,which(!colnames(seqTab) %in% c("decipherCluster","Row.names.y"))],
                     list(seqTab$decipherCluster),function(x){
                        if(class(x) %in% c("integer","numeric")){
                            out <- sum(x)
                        }else{
                            out <- paste(unique(x),sep="|",collapse="|")
                        }
                        out
                        })
colnames(clusTab)[1:2] <- c("decipherCluster","OTUs")


print("writing output")
write.table(clusTab,
            snakemake@output[[2]],sep="\t",col.names=NA,quote=F)
saveRDS(clusTab,
        snakemake@output[[3]])

