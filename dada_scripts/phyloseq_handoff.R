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
library(phyloseq)
library(Biostrings)
library(ape)

# File parsing
seqTab <- readRDS(snakemake@input[[1]])
sInfo <- read.delim(snakemake@input[[2]],stringsAsFactors=F,row.names=1)

 seqMat <- seqTab[,colnames(seqTab) %in% rownames(sInfo)]
 rownames(seqMat) <- seqTab$Row.names
 taxMat <- as.matrix(seqTab[,!colnames(seqTab) %in% c(rownames(sInfo),"Row.names")])
 row.names(taxMat) <- seqTab$Row.names
 seqPhy <- phyloseq(otu_table(seqMat,taxa_are_rows = T),
                sample_data(sInfo),
                tax_table(as.matrix(taxMat)))

 seqs <- Biostrings::DNAStringSet(seqTab$Row.names)
 names(seqs) <- seqTab$Row.names
 seqPhy <- merge_phyloseq(seqPhy, seqs)
 taxa_names(seqPhy) <- tax_table(seqPhy)[,"OTU"]
if(snakemake@params[["currentStep"]]=="post"&snakemake@config[['postprocessing']][['treeing']]){
 tree <- read.tree(snakemake@input[[3]])
 seqPhy <- phyloseq(otu_table(seqPhy),
                    sample_data(seqPhy),
                    tax_table(seqPhy),
                    phy_tree(tree))
}
writeRDS(seqPhy,snakemake@output[[1]])


