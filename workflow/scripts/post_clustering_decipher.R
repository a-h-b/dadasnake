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
if(grepl("OTU_",names(seqs)[1])){
  names(seqs) <- gsub("OTU_","ASV_",names(seqs))
}

print(paste0("loading ASV table",snakemake@input[[2]]))
seqTab <- readRDS(snakemake@input[[2]])
if(any(colnames(seqTab)=="OTU")){
 seqTab$OTU <- gsub("OTU_","ASV_",seqTab$OTU)
 colnames(seqTab)[colnames(seqTab)=="OTU"] <- "ASV"
}

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
colnames(clusters) <- "OTU"


write.table(clusters,
            snakemake@output[[1]],sep="\t",col.names=NA,quote=F)

print("calculate clustered table")
seqTab <- merge(clusters,seqTab,by.x=0,by.y="ASV")

clusSeqs <- aggregate(seqTab$Row.names,
                      list(seqTab$OTU),function(x){
                        x[1]
                      })
cluseq <- seqs[clusSeqs$x]
num_digits <- as.character(ceiling(log10(nrow(clusSeqs))))
clusSeqs$Group.1 <- sprintf(paste0("OTU_%0",num_digits,"d"),clusSeqs$Group.1)
names(cluseq) <- clusSeqs$Group.1

writeXStringSet(cluseq,snakemake@output[[2]])

clusTab <- aggregate(seqTab[,which(!colnames(seqTab) %in% c("OTU","Row.names.y"))],
                     list(seqTab$OTU),function(x){
                        if(class(x) %in% c("integer","numeric")){
                            out <- sum(x)
                        }else{
                            out <- paste(unique(x),sep="|",collapse="|")
                        }
                        out
                        })
colnames(clusTab)[which(sapply(clusTab,class)=="character")] <- paste0("consensus.",colnames(clusTab)[which(sapply(clusTab,class)=="character")]) 
colnames(clusTab)[1:2] <- c("OTU","ASV")
clusTab$OTU <- sprintf(paste0("OTU_%0",num_digits,"d"),clusTab$OTU)
clusTab$Row.names <- sapply(clusTab$OTU, function(x) as.character(cluseq[[x]]))

print("writing output")
write.table(clusTab,
            snakemake@output[[3]],sep="\t",col.names=NA,quote=F)
saveRDS(clusTab,
        snakemake@output[[4]])

