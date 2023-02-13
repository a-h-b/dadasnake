log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(Biostrings)

print(paste0("loading ASV table",snakemake@input[[2]]))
seqTab <- readRDS(snakemake@input[[2]])
if(any(colnames(seqTab)=="OTU")){
 seqTab$OTU <- gsub("OTU_","ASV_",seqTab$OTU)
 colnames(seqTab)[colnames(seqTab)=="OTU"] <- "ASV"
}

print("loading clusters")
clusters <- read.delim(snakemake@input[[3]],sep="\t", header=F)[,c(1,2,9)]
clusters <- clusters[clusters$V1 != "S",c(2,3)]
colnames(clusters) <- c("OTU","ASV")

print(paste0("loading cluster sequence ",snakemake@input[[1]]))
clusSeqs <- readDNAStringSet(snakemake@input[[1]])

num_digits <- as.character(max(ceiling(log10(length(unique(clusters$OTU)))),nchar(names(clusSeqs)[1])-4))
clusters$OTU <- sprintf(paste0("OTU_%0",num_digits,"d"),clusters$OTU)


print("calculate clustered table")
seqTab <- merge(clusters,seqTab,by="ASV")

print(dim(seqTab))

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

