log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(Biostrings)

print(paste0("loading ASV table",snakemake@input[[1]]))
seqTab <- readRDS(snakemake@input[[1]])

print("loading clusters")
clusters <- read.delim(snakemake@input[[2]],sep="\t", header=F,comment.char="S")[,c(2,9)]
colnames(clusters) <- c("vsearchCluster","OTU")

print("calculate clustered table")
seqTab <- merge(clusters,seqTab,by="OTU")

clusTab <- aggregate(seqTab[,which(!colnames(seqTab) %in% c("vsearchCluster","Row.names"))],
                     list(seqTab$vsearchCluster),function(x){
                        if(class(x) %in% c("integer","numeric")){
                            out <- sum(x)
                        }else{
                            out <- paste(unique(x),sep="|",collapse="|")
                        }
                        out
                        })
colnames(clusTab)[1:2] <- c("vsearchCluster","OTUs")


print("writing output")
write.table(clusTab,
            snakemake@output[[1]],sep="\t",col.names=NA,quote=F)
saveRDS(clusTab,
        snakemake@output[[2]])

