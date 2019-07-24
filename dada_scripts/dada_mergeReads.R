args<-commandArgs(TRUE)
option <- as.numeric(args[1]) #which parameter set to use
outfile <- args[2]

library(BiocParallel)
register(SerialParam())
library(dada2)
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")

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
  ggplot(data = plotdf, aes(x = Cycle, y = Score)) + geom_tile(aes(fill = Count)) + 
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
                                                                                              strip.text.x = element_blank())
}


# File parsing
path <- getwd() # CHANGE ME to the directory containing your demultiplexed forward-read fastq files
fastqFs <- sort(list.files(paste0(path,"/cutraw"), pattern="fwd.fastq"))
fastqRs <- sort(list.files(paste0(path,"/cutraw"), pattern="rvs.fastq"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

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

save.image(paste0(filtpath,"/filterWS.Rdata"))
#print("saving filter stats")
#saveRDS(fileout,paste0(filtpath,"/filterstats.RDS"))
pdf(paste0(filtpath,"/QC_filt.pdf"),
    width=8,height=11,pointsize=7)
plotQualityProfile(filtFs)
plotQualityProfile(filtRs)
dev.off()


sample.namesR <- gsub(".rvs.fastq","",fastqRs) 
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

print("learning error models F")
# Learn forward error rates
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=F)
errF <- dadaFs.lrn[[1]]$err_out
pdf(paste0(filtpath,"/errorModels_fwd.pdf"),
    width=8,height=11,pointsize=7)
plotErrors(dadaFs.lrn[[1]], nominalQ=TRUE)
dev.off()
saveRDS(errF,paste0(filtpath,"/errorModels_fwd.RDS"))
rm(errF)
print("learning error models R")
# Learn reverse error rates
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=F)
errR <- dadaRs.lrn[[1]]$err_out
pdf(paste0(filtpath,"/errorModels_rvs.pdf"),
    width=8,height=11,pointsize=7)
plotErrors(dadaRs.lrn[[1]], nominalQ=TRUE)
dev.off()
saveRDS(errR,paste0(filtpath,"/errorModels_rvs.RDS"))

print("merging")
errF <- readRDS(paste0(filtpath,"/errorModels_fwd.RDS"))
errR <- readRDS(paste0(filtpath,"/errorModels_rvs.RDS"))
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=FALSE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=FALSE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
print("saving sequence tab")
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, paste0(filtpath,"/seqtab.rds"))
write.table(paste0(filtpath,"/seqtab.rds"),outfile,append=T,row.names=F,col.names=F,quote=F)
