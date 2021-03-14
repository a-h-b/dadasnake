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
library(phyloseq)
library(Biostrings)
library(ape)

# File parsing
if(file.info(snakemake@input[[1]])$size == 0){
  print("no reads in OTU table.")
  system(paste0("touch ",snakemake@output[[1]]))
}else{
  seqTab <- readRDS(snakemake@input[[1]])
  sInfo <- read.delim(snakemake@input[[2]],stringsAsFactors=F,row.names=1)
  rownames(sInfo) <- gsub("-",".",rownames(sInfo))
  if(nrow(sInfo)>1){
    seqMat <- as.matrix(seqTab[,colnames(seqTab) %in% ifelse(grepl("^[[:digit:]]",rownames(sInfo)),
                                                             paste0("X",rownames(sInfo)),
                                                              rownames(sInfo))])
    rownames(seqMat) <- seqTab$Row.names
    if(any(!colnames(seqTab) %in% c(ifelse(grepl("^[[:digit:]]",rownames(sInfo)),
                                                             paste0("X",rownames(sInfo)),
                                                              rownames(sInfo)),"Row.names"))){
      taxMat <- as.matrix(seqTab[,!colnames(seqTab) %in% c(ifelse(grepl("^[[:digit:]]",rownames(sInfo)),
                                                             paste0("X",rownames(sInfo)),
                                                              rownames(sInfo)),"Row.names")])
    }else{
      taxMat <- matrix(seqTab$Row.names,ncol=1,nrow=length(seqTab$Row.names))
    }
    row.names(taxMat) <- seqTab$Row.names
    if(is.null(colnames(taxMat))) colnames(taxMat) <- "OTU"
    if(any(grepl("^[[:digit:]]",rownames(sInfo)))) rownames(sInfo) <- ifelse(grepl("^[[:digit:]]",rownames(sInfo)),
                                                             paste0("X",rownames(sInfo)),
                                                              rownames(sInfo))
    sInfo <- sInfo[ifelse(grepl("^[[:digit:]]",rownames(sInfo)),
                         paste0("X",rownames(sInfo)),
                         rownames(sInfo)) %in% colnames(seqMat),]
    seqPhy <- phyloseq(otu_table(seqMat,taxa_are_rows = T),
                 sample_data(sInfo),
                 tax_table(taxMat))

    seqs <- Biostrings::DNAStringSet(seqTab$Row.names)
    names(seqs) <- seqTab$Row.names
    seqPhy <- merge_phyloseq(seqPhy, seqs)
    taxa_names(seqPhy) <- tax_table(seqPhy)[,"OTU"]
    if(snakemake@params[["currentStep"]]=="post"&snakemake@config[['postprocessing']][['treeing']][['do']]){
       tree <- read.tree(snakemake@input[[3]])
       seqPhy <- phyloseq(otu_table(seqPhy),
                      sample_data(seqPhy),
                      tax_table(seqPhy),
                      phy_tree(tree))
    }
    saveRDS(seqPhy,snakemake@output[[1]])
  }else{
    saveRDS(NULL,snakemake@output[[1]])
  }
}
