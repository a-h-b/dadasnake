args<-commandArgs(TRUE)
#option <- as.numeric(args[1]) #which parameter set to use
infile <- args[1]
outpath <- args[2]

library(BiocParallel)
register(SerialParam())
library(dada2)
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(Biostrings)

#new mergeSequenceTable function
mergeSequenceTables <- function(tables, orderBy = "abundance"){
    sample.names <- rownames(tables[[1]])
    for (i in seq(2, length(tables))) {
        sample.names <- c(sample.names, rownames(tables[[i]]))
    }
    if (any(duplicated(sample.names))) {
        stop("Duplicated sample names detected in the rownames.")
    }
    seqs <- unique(c(sapply(tables, colnames), recursive = TRUE))
    rval <- matrix(0L, nrow = length(sample.names), ncol = length(seqs))
    rownames(rval) <- sample.names
    colnames(rval) <- seqs
    for (tab in tables) {
        rval[rownames(tab), colnames(tab)] <- tab
    }
    if (!is.null(orderBy)) {
        if (orderBy == "abundance") {
            rval <- rval[, order(colSums(rval), decreasing = TRUE), 
                drop = FALSE]
        }
        else if (orderBy == "nsamples") {
            rval <- rval[, order(colSums(rval > 0), decreasing = TRUE), 
                drop = FALSE]
        }
    }
    rval
}


path <- getwd()


print("merging runs")
# # Merge multiple runs (if necessary)
seqTabList <- read.delim(infile,header=F,stringsAsFactors=F)[,1]
if(length(seqTabList)>1){
seqTabs <- list()
for(i in 1:length(seqTabList)) seqTabs[[i]] <- readRDS(seqTabList[i])
seqtab <- mergeSequenceTables(tables=seqTabs)
}else{
seqtab <- readRDS(seqTabList[1])
}
print("Removing chimeras")
seqtab.nochime <- removeBimeraDenovo(seqtab, tableMethod="consensus")
print("saving")
#saveRDS(seqtab.nochime, paste0(outpath,"/seqtab.nochime.rds"))
# Assign taxonomy
#seqtab.nochime2 <- seqtab.nochime
#seqs <- reverseComplement(DNAStringSet(colnames(seqtab.nochime)))
seqs <- DNAStringSet(colnames(seqtab.nochime))
names(seqs) <- sprintf("OTU_%06d",1:length(seqs))
writeXStringSet(seqs,paste0(outpath,"/seqs.nochime.fasta"))


