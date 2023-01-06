log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(Biostrings)

print(paste0("loading ASV table",snakemake@input[[1]]))
seqTab <- readRDS(snakemake@input[[1]])
if(any(colnames(seqTab)=="OTU")){
 seqTab$OTU <- gsub("OTU_","ASV_",seqTab$OTU)
 colnames(seqTab)[colnames(seqTab)=="OTU"] <- "ASV"
}

print("loading clusters")
clusters <- read.delim(snakemake@input[[2]],sep="\t", header=F,comment.char="S")[,c(2,9)]
colnames(clusters) <- c("OTU","ASV")

print(paste0("loading cluster sequence ",snakemake@input[[3]]))
clusSeqs <- readDNAStringSet(snakemake@input[[3]])


print("calculate clustered table")
seqTab <- merge(clusters,seqTab,by="ASV")

clusTab <- aggregate(seqTab[,which(!colnames(seqTab) %in% c("OTU","Row.names"))],
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
clusTab$Row.names <- sapply(clusTab$OTU, function(x) as.character(clusSeqs[[x]]))

print("writing output")
write.table(clusTab,
            snakemake@output[[1]],sep="\t",col.names=NA,quote=F)
saveRDS(clusTab,
        snakemake@output[[2]])

