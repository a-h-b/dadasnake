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
library(dada2)


# File parsing
fastqF <- snakemake@input[[1]]
fastqR <- snakemake@input[[2]]
filtF <- snakemake@output[[1]]
filtR <- snakemake@output[[2]]
sampleName <- gsub("preprocessing/","run.",gsub("/","",gsub(".fwd.fastq","",fastqF)))
filtpath <- gsub("/.+","",filtF)

#output table with number of reads pre-filtering

# shorter/longer truncation; no truncation, truncQ higher; maxEE = 1
partab <- data.frame("option"=1:20,
                     "filtdir"=paste0(path,"/",c("option."),sprintf("%02d",1:20)),
                     "maxEE"=c(rep(2,3),rep(2,5),rep(c(1,0.7),each=3),rep(c(1,0.7),each=3)),
                      "truncQ"=c(rep(10,3),25,20,15,10,2,rep(c(25,10,2),times=2),rep(c(25,10,2),times=2)),
                     "len1"=c(90,100,130,rep(0,11),rep(100,6)),
                     "len2"=c(80,90,120,rep(0,11),rep(90,6)),
                     stringsAsFactors=F)

filtpath <- file.path(partab$filtdir[partab$option==option]) # Filtered files go into this subdirectory

print("filtering")

if(!file_test("-d", filtpath)) dir.create(filtpath)
sample.names <- gsub(".fwd.fastq","",fastqFs)
filtFs <- file.path(filtpath, paste0(sample.names, "_fwd_filt.fastq.gz"))
filtRs <- file.path(filtpath, paste0(sample.names, "_rvs_filt.fastq.gz"))
for(i in seq_along(fastqFs)) {
  fastqPairedFilter(c(paste0("cutraw/",fastqFs[i]),paste0("cutraw/",fastqRs[i])), 
                    c(filtFs[i],filtRs[i]),
                    truncLen=c(partab$len1[partab$option==option],partab$len2[partab$option==option]), 
                    maxEE=partab$maxEE[partab$option==option], 
                    truncQ=partab$truncQ[partab$option==option], maxN=0, rm.phix=TRUE,
                        compress=TRUE, verbose=TRUE)
}

#output table with number of reads post-filtering
print("filtering done")
