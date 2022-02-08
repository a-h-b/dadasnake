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


sampleFile <- snakemake@input[[1]]
filts <- unique(unlist(snakemake@input)[-1])

print(paste0("reading sample table from ",sampleFile))
sampleTab <- read.delim(sampleFile,stringsAsFactors=F)

crun <- unlist(strsplit(filts[1],split="/"))[2]
csam <- gsub(".fwd.fastq.gz","", unlist(strsplit(filts[1],split="/"))[3])

creads <- sampleTab$reads_filtered_fwd[sampleTab$sample==csam&sampleTab$run==crun][1]

if(!snakemake@config[['downsampling']][['use_total']]){
  downsize <- as.numeric(snakemake@config[['downsampling']][['number']])
  if(is.na(downsize)) downsize <- 0 else downsize <- as.numeric(downsize)
  if(creads > 0){
    sreads <- sapply(unique(sampleTab$sample),
                 function(y) sum(sapply(unique(sampleTab$run[sampleTab$sample==y]),
                 function(x) sampleTab$reads_filtered_fwd[sampleTab$sample==y&sampleTab$run==x][1])))
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
  }else{
     print(paste0("no reads in ",crun,"/",csam))
     system(paste("touch", snakemake@output[[1]]))
     system(paste("touch", snakemake@output[[2]]))
  }
}else{
  all_down <- as.numeric(snakemake@config[['downsampling']][['total_number']])
  if(all_down<=0){
    print("Invalid total number of reads, no downsampling was performed")
    system(paste("cp ",filts[1]," ",snakemake@output[[1]]))
    system(paste("cp ",filts[2]," ",snakemake@output[[2]]))
  }else{
    all_run <- sum(sampleTab$reads_filtered_fwd[sampleTab$run==crun])
    downsize <- ceiling(creads * all_run/all_down)
  }
}

if(downsize>0){
    system(paste("seqtk sample -s",snakemake@config[['downsampling']][['seed']],filts[1],downsize," | gzip >",snakemake@output[[1]]))
    system(paste("seqtk sample -s",snakemake@config[['downsampling']][['seed']],filts[2],downsize," | gzip >",snakemake@output[[2]]))
}else{
    print(paste0("not enough reads in ",csam))
    system(paste("touch", snakemake@output[[1]]))
    system(paste("touch", snakemake@output[[2]]))
}

print("done")
