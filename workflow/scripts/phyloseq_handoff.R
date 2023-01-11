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
  if(snakemake@params[["what"]]=="ASV"){
    units <- "ASV"
  }else{
    units <- "OTU"
  }
  seqTab <- readRDS(snakemake@input[[1]])
  if(units=="ASV" & any(colnames(seqTab)=="OTU")){
    seqTab$OTU <- gsub("OTU_","ASV_",seqTab$OTU)
    colnames(seqTab)[which(colnames(seqTab)=="OTU")] <- "ASV"
  }
  sInfo <- read.delim(snakemake@input[[2]],stringsAsFactors=F,row.names=1)
  if(!any(grepl("-",colnames(seqTab)))) rownames(sInfo) <- gsub("-",".",rownames(sInfo))
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

    if(is.null(colnames(taxMat))){
      colnames(taxMat) <- units
    }
    if(snakemake@params[['cluster']]=="add" & snakemake@params[['what']]=="ASV"){
      cFile <- which(sapply(snakemake@input,function(x) grepl("cluster",x)))
      if(snakemake@params[['clusterMeth']]=="vsearch"){
        clusterInfo <- read.delim(snakemake@input[[cFile]],sep="\t", header=F)[,c(1,2,9)]
        clusterInfo <- clusterInfo[clusterInfo$V1 != "S",c(2,3)]
        colnames(clusterInfo) <- c("OTU","ASV")
        num_digits <- as.character(ceiling(log10(length(unique(clusterInfo$OTU)))))
        clusterInfo$OTU <- sprintf(paste0("OTU_%0",num_digits,"d"),clusterInfo$OTU)
      }else{
        clusterInfo <- read.delim(snakemake@input[[cFile]],sep="\t")
        colnames(clusterInfo) <- c("ASV","OTU")
      }

      taxMat <- cbind(taxMat, sapply(taxMat[,which(colnames(taxMat)=="ASV")],function(x) clusterInfo$OTU[clusterInfo$ASV==x]))
      colnames(taxMat)[ncol(taxMat)] <- "OTU"
    }
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
    taxa_names(seqPhy) <- tax_table(seqPhy)[,units]
    if(snakemake@params[['tree']]=="add"){
       tree <- read.tree(snakemake@input[[3]])
       if(snakemake@params[['what']]=="ASV" & grepl("OTU_",tree$tip.label)) tree$tip.label <- gsub("OTU_","ASV_", tree$tip.label)
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
