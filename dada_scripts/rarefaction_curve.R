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
library(vegan)

### rarefaction curve function
dadasnake_rarecurve <- function (x, step = 100, sample, xlab = "reads", ylab = "sequence variants", 
                                 label = TRUE, col="black", ...) {
  tot <- rowSums(x)
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
  if(length(col)<nrow(x)) col <- rep(col,nrow(x))
  par(mar=c(3.3,4.5,0.3,0.2),mgp=c(3.5,0.6,0))
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "", ylab = ylab, 
       type = "n", las=1, ...)
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
print("reading input")
seqTab <- readRDS(snakemake@input[[1]])
sInfo <- read.delim(snakemake@input[[2]],stringsAsFactors=F,row.names=1)

print("formatting OTU table")
seqMat <- as.matrix(seqTab[,colnames(seqTab) %in% rownames(sInfo)])
rownames(seqMat) <- seqTab$Row.names

if(any(colSums(seqMat)<1)){
  seqMat <- seqMat[,colSums(seqMat)>0]
  print("Some samples contained no reads and are omitted from rarefaction curve:")
  print(colnames(seqMat)[colSums(seqMat)<1])
}

print("plotting rarefaction curves")
pdf(snakemake@output[[1]],width=15/2.54,height = 12/2.54,pointsize = 7)
dadasnake_rarecurve(t(seqMat),label=F,lty=1)
dev.off()

