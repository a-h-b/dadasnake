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
if(!require(dada2)){
  BiocManager::install("GenomeInfoDbData",update=F,ask=F)
  require(dada2)
}
library(ShortRead)
library(ggplot2)

#new QC function
plotQualityProfile <- function(fl, n = 5e+05){
  statdf <- data.frame(Cycle = integer(0), Mean = numeric(0), 
                       Q25 = numeric(0), Q50 = numeric(0), Q75 = numeric(0), 
                       file = character(0))
  anndf <- data.frame(minScore = numeric(0), label = character(0), 
                      rclabel = character(0), file = character(0))
  FIRST <- TRUE
  for (f in fl) {
    srqa <- qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    rc <- srqa[["readCounts"]]$read
    if (rc >= n) {
      rclabel <- paste("Reads >= ", n)
    }
    else {
      rclabel <- paste("Reads: ", rc)
    }
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count, 
                                                          df$Cycle)
    get_quant <- function(xx, yy, q) {
      xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
    }
    q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.25), simplify = TRUE)
    q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.5), simplify = TRUE)
    q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.75), simplify = TRUE)
    if (!all(sapply(list(names(q25s), names(q50s), names(q75s)), 
                    identical, rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    if (FIRST) {
      plotdf <- cbind(df, file = f)
      FIRST <- FALSE
    }
    else {
      plotdf <- rbind(plotdf, cbind(df, file = f))
    }
    statdf <- rbind(statdf, data.frame(Cycle = as.integer(rownames(means)), 
                                       Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s), 
                                       Q75 = as.vector(q75s), file = f))
    anndf <- rbind(anndf, data.frame(minScore = min(df$Score), 
                                     label = basename(f), rclabel = rclabel, file = f))
  }
  anndf$minScore <- min(anndf$minScore)
  print(ggplot(data = plotdf, aes(x = Cycle, y = Score)) + geom_tile(aes(fill = Count)) + 
    scale_fill_gradient(low = "#F5F5F5", high = "black") + 
    geom_line(data = statdf, aes(y = Mean), color = "#66C2A5") + 
    geom_line(data = statdf, aes(y = Q25), color = "#FC8D62", 
              size = 0.25, linetype = "dashed") + geom_line(data = statdf, 
                                                            aes(y = Q50), color = "#FC8D62", size = 0.25) + geom_line(data = statdf, 
                                                                                                                      aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
    ylab("Quality Score") + xlab("Cycle") + theme_bw() + 
    theme(panel.grid = element_blank()) + guides(fill = FALSE) + 
    geom_text(data = anndf, aes(x = minScore + 2, label = label), 
              y = 1, hjust = 0, vjust = 0) + geom_text(data = anndf, 
                                                       aes(x = minScore + 2, label = rclabel), y = 1, hjust = 0, 
                                                       vjust = 2) + facet_wrap(~file) + theme(strip.background = element_blank(), 
                                                                                              strip.text.x = element_blank()))
}


# File parsing
path <- snakemake@params[['path']] 
fastqFs <- sort(list.files(path, pattern="fwd.fastq"))
fastqRs <- sort(list.files(path, pattern="rvs.fastq"))
sizeFs <- sapply(paste0(path,"/",fastqFs),function(x){
            if(grepl(".gz$",x)){
             as.numeric(unlist(strsplit(system2("zcat",args=c(x," 2>/dev/null | head | wc -l"),stdout=T),split=" "))[1])
            }else{
             file.info(x)$size
            }
          })
fastqFs <- fastqFs[sizeFs>0]
sizeRs <- sapply(paste0(path,"/",fastqRs),function(x) {
            if(grepl(".gz$",x)){
             as.numeric(unlist(strsplit(system2("zcat",args=c(x," 2>/dev/null | head | wc -l"),stdout=T),split=" "))[1])
            }else{
             file.info(x)$size
            }
          })
fastqRs <- fastqRs[sizeRs>0]

# QC
pdf(snakemake@output[[1]],
    width=8,height=11,pointsize=7)
if(length(fastqFs)>0){
  if(floor(length(fastqFs)/36)>0){
    for(i in 1:floor(length(fastqFs)/36)){
      plotQualityProfile(paste0(path,"/",fastqFs[(36*(i-1))+1:36]))
    }
  }
  if(length(fastqFs) %% 36 > 0){
    plotQualityProfile(paste0(path,"/",fastqFs[(36*(floor(length(fastqFs)/36))+1):length(fastqFs)]))
  }
}
dev.off()

pdf(snakemake@output[[2]],
    width=8,height=11,pointsize=7)
if(length(fastqRs)>0){
  if(floor(length(fastqRs)/36)>0){
    for(i in 1:floor(length(fastqRs)/36)){
      plotQualityProfile(paste0(path,"/",fastqRs[(36*(i-1))+1:36]))
    }
  }
  if(length(fastqRs) %% 36 > 0){
    plotQualityProfile(paste0(path,"/",fastqRs[(36*(floor(length(fastqRs)/36))+1):length(fastqRs)]))
  }
}
dev.off()


