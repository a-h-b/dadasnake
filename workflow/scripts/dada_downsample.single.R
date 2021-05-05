log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(BiocParallel)
if (snakemake@threads > 1) {
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}else{
    parallel <- FALSE
    register(SerialParam())
}
if(!require(ShortRead)){
  BiocManager::install("GenomeInfoDbData",update=F,ask=F)
  library(ShortRead)
}


sampleFile <- snakemake@input[[1]]
filts <- unique(unlist(snakemake@input)[-1])

print(paste0("reading sample table from ",sampleFile))
sampleTab <- read.delim(sampleFile,stringsAsFactors=F)

downsize <- snakemake@config[['downsampling']][['number']]
if(is.na(as.numeric(downsize))) downsize <- 0 else downsize <- as.numeric(downsize)

crun <- unlist(strsplit(filts[1],split="/"))[2]
csam <- gsub(".fastq.gz","", unlist(strsplit(filts[1],split="/"))[3])

creads <- sampleTab$reads_filtered[sampleTab$sample==csam&sampleTab$run==crun][1]

if(creads > 0){
  sreads <- sapply(unique(sampleTab$sample),
                 function(y) sum(sapply(unique(sampleTab$run[sampleTab$sample==y]),
                 function(x) sampleTab$reads_filtered[sampleTab$sample==y&sampleTab$run==x][1])))
  treads <- sreads[csam]
  if(snakemake@config[['downsampling']][['min']]) downsize <- max(downsize,min(sreads[sreads>0]))
  if(treads >= downsize){
    if(length(unique(sampleTab$run[sampleTab$sample==csam]))>1){
      cpart <- creads/treads
      downsize <- round(downsize*cpart) 
    } 
  }else{
    downsize <- 0
  } 
  if(downsize>0){
    set.seed(snakemake@config[['downsampling']][['seed']])
  
    sampler <- FastqSampler(filts[1],downsize)
  
    writeFastq(yield(sampler),snakemake@output[[1]],compress=T)
  }else{
    print(paste0("not enough reads in ",csam))
    system(paste("touch", snakemake@output[[1]]))
  }  

}else{
 print(paste0("no reads in ",crun,"/",csam))
 system(paste("touch", snakemake@output[[1]]))
}

print("done")
