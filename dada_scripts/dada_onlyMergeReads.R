log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

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
library(dada2)
#library(ggplot2)

# File parsing
errFfile <- snakemake@input[[1]]
errRfile <- snakemake@input[[2]]

#filtFs <- sort(grep("fwd.fastq",unlist(snakemake@input),value=T))
#filtRs <- sort(grep("rvs.fastq",unlist(snakemake@input),value=T))
filtF <- snakemake@input[[3]]
filtR <- snakemake@input[[4]]

#if(length(filtFs) != length(filtRs)) stop("Forward and reverse files do not match.")

sampleNameF <- gsub("/","",gsub("filtered/","run.",gsub(".fwd.fastq.gz","",filtF)))
sampleNameR <- gsub("/","",gsub("filtered/","run.",gsub(".rvs.fastq.gz","",filtR)))

if(sampleNameF!=sampleNameR) stop("Forward and reverse files do not match.")

sampleName <- gsub(".+/","",gsub(".fwd.fastq.gz","",filtF))

names(filtF) <- sampleName
names(filtR) <- sampleName

#dadafile <- snakemake@output[[1]]
#dadatab <- snakemake@output[[2]]
#dadatabtsv <- snakemake@output[[3]]
mergefile <- snakemake@output[[1]]

print("merging")
errF <- readRDS(errFfile)
errR <- readRDS(errRfile)

# Sample inference and merger of paired-end reads
print(paste0("make dada object and merge, ",sampleName))
#mergers <- vector("list", length(sampleNames))
#names(mergers) <- sampleNames
#for(sam in sampleNames) {
#  print(paste0("Processing: ", sam))
  derepF <- derepFastq(filtF)
  dadaF <- dada(derepF, err=errF, multithread=snakemake@threads)
  derepR <- derepFastq(filtR)
  dadaR <- dada(derepR, err=errR, multithread=snakemake@threads)
  merger <- mergePairs(dadaF, derepF, dadaR, derepR,
            minOverlap=snakemake@config[['pair_merging']][['min_overlap']],
            maxMismatch=snakemake@config[['pair_merging']][['max_mismatch']],
            justConcatenate=snakemake@config[['pair_merging']][['just_concatenate']],
            trimOverhang=snakemake@config[['pair_merging']][['trim_overhang']])
#  mergers[[sam]] <- merger
#}
saveRDS(merger,mergefile)
#print("making sequence tab")
# Construct sequence table
#seqtab <- makeSequenceTable(mergers)
#saveRDS(seqtab, dadatab)
#write.table(t(seqtab),dadatabtsv,col.names=NA,sep="\t",quote=F)

print("done")
