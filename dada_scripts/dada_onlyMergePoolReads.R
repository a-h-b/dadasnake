log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(BiocParallel)
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}else{
    parallel <- FALSE
    register(SerialParam())
}
library(dada2)

# File parsing
dadaFfile <- snakemake@input[[1]]
dadaRfile <- snakemake@input[[2]]
derepFfile <- snakemake@input[[3]]
derepRfile <- snakemake@input[[4]]
mergefile <- snakemake@output[[1]]

dadaF <- readRDS(dadaFfile)
dadaR <- readRDS(dadaRfile)

if(names(dadaF)!=names(dadaR)) stop("Forward and reverse files do not match.")

derepF <- readRDS(derepFfile)
derepR <- readRDS(derepRfile)

print("merging")
# merger of paired-end reads
print("merge reads")
merger <- mergePairs(dadaF, derepF, dadaR, derepR,
                     minOverlap=snakemake@config[['pair_merging']][['min_overlap']],
                     maxMismatch=snakemake@config[['pair_merging']][['max_mismatch']],
                     justConcatenate=snakemake@config[['pair_merging']][['just_concatenate']],
                     trimOverhang=snakemake@config[['pair_merging']][['trim_overhang']])
saveRDS(merger,mergefile)

print("done")
