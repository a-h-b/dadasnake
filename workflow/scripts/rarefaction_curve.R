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
library(vegan)

### rarefaction curve function
dadasnake_rarecurve <- function (x, step = 100, sample, xlab = "reads", ylab = snakemake@params[["thing"]], 
                                 ymax,xmax,
                                 label = TRUE, col="black", ...) {
  tot <- rowSums(x)
  if(any(tot < step)){
    step <- min(tot)
    print(paste0("New step size: ",step))
  }
  S <- specnumber(x)
  nr <- nrow(x)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  if(missing(ymax)) ymax <- max(Smax)
  if(missing(xmax)) xmax <- max(Nmax)
  if(length(col)<nrow(x)) col <- rep(col,nrow(x))
  par(mar=c(3.3,4.5,0.3,0.2),mgp=c(3.5,0.6,0))
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "", ylab = ylab, 
       type = "n", las=1, ylim=c(1,ymax), xlim=c(1,xmax), ...)
  mtext(xlab,1,2)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_len(length(out))) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col=col[ln],...)
  }
  if (label) {
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}


# File parsing
if(file.info(snakemake@input[[1]])$size == 0){
  print("no reads in OTU table.")
  system(paste0("touch ",snakemake@output[[1]]))
}else{
  print("reading input")
  seqTab <- readRDS(snakemake@input[[1]])
  sInfo <- read.delim(snakemake@input[[2]],stringsAsFactors=F,row.names=1)
  if((ncol(seqTab)==3 &!colnames(seqTab)[2] %in% rownames(sInfo)) | 
       !colnames(seqTab)[3] %in% ifelse(grepl("^[[:digit:]]",c(rownames(sInfo),sInfo$sample)),
                                        paste0("X",c(rownames(sInfo),sInfo$sample)),
                                        c(rownames(sInfo),sInfo$sample))) colnames(seqTab)[which(colnames(seqTab)=="V1")] <- rownames(sInfo)

  sams <- which(colnames(seqTab) %in% ifelse(grepl("^[[:digit:]]",c(rownames(sInfo),sInfo$sample)),
                                        paste0("X",c(rownames(sInfo),sInfo$sample)), 
                                        c(rownames(sInfo),sInfo$sample)))
  samno <- length(sams)
  if(samno>0){
    if(samno > 1){
      ymax <- max(apply(seqTab[,sams],2,function(x) length(which(x>0))))
      xmax <- max(colSums(seqTab[,sams]))
    }else{
      ymax <- length(which(seqTab[,sams]>0))
      xmax <- sum(seqTab[,sams])
    }
    pdf(snakemake@output[[1]],width=15/2.54,height = 12/2.54,pointsize = 7)
    if(floor(samno/72)>0){
      for(i in 1:floor(samno/72)){
        print(i) 
        print("formatting OTU table")
        seqMat <- as.matrix(seqTab[,sams[(72*(i-1))+1:72]])
        rownames(seqMat) <- seqTab$Row.names

        if(any(colSums(seqMat)<1)){
          print("Some samples contained no reads and are omitted from rarefaction curve:")
          print(colnames(seqMat)[colSums(seqMat)<1])
          seqMat <- seqMat[,colSums(seqMat)>0]
        }

        print("plotting rarefaction curves")
        dadasnake_rarecurve(t(seqMat),label=F,lty=1,ymax=ymax,xmax=xmax)
      }
    }else{
     i <- 0
    }
    if(samno %% 72 > 0){
      print(i+1)
      seqMat <- as.matrix(seqTab[,sams[(72*(floor(length(sams)/72))+1):length(sams)]])
      rownames(seqMat) <- seqTab$Row.names
      if(any(colSums(seqMat)<1)){
          print("Some samples contained no reads and are omitted from rarefaction curve:")
          print(colnames(seqMat)[colSums(seqMat)<1])
          seqMat <- seqMat[,colSums(seqMat)>0]
        }

        print("plotting rarefaction curves")
      dadasnake_rarecurve(t(seqMat),label=F,lty=1,ymax=ymax,xmax=xmax)
    }
    dev.off()
  }else{
    print("No columns to rarefy")
    system(paste0("touch ",snakemake@output[[1]]))
  }
}
