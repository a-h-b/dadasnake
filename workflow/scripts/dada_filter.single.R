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

fastqFilter0 <- function (fn, fout, truncQ = 2, truncLen = 0, maxLen = Inf, 
          minLen = 20, trimLeft = 0, trimRight = 0, maxN = 0, minQ = 0, 
          maxEE = Inf, rm.phix = TRUE, rm.lowcomplex = 0, orient.fwd = NULL, 
          n = 1e+06, OMP = TRUE, qualityType = "Auto", compress = TRUE, 
          verbose = FALSE, ...) {
  if (!OMP) {
    ompthreads <- .Call(ShortRead:::.set_omp_threads, 1L)
    on.exit(.Call(ShortRead:::.set_omp_threads, ompthreads))
  }
  if (any(sapply(list(truncQ, truncLen, maxLen, minLen, trimLeft, 
                      trimRight, maxN, minQ, maxEE), length) > 1)) {
    stop("Filtering and trimming arguments should be of length 1 when processing single-end (rather than paired-end) data.")
  }
  start <- max(1, trimLeft + 1, na.rm = TRUE)
  end <- truncLen
  if (end < start) {
    end = NA
  }
  end <- end - start + 1
  if (fn == fout) {
    stop("The output and input files must be different.")
  }
  f <- FastqStreamer(fn, n = n)
  on.exit(close(f))
  if (file.exists(fout)) {
    if (file.remove(fout)) {
      if (verbose) 
        message("Overwriting file:", fout)
    }
    else {
      stop("Failed to overwrite file:", fout)
    }
  }
  first = TRUE
  inseqs = 0
  outseqs = 0
  while (length(suppressWarnings(fq <- yield(f, qualityType = qualityType)))) {
    inseqs <- inseqs + length(fq)
    if (!is.null(orient.fwd)) {
      if (!C_isACGT(orient.fwd)) 
        stop("Non-ACGT characters detected in orient.fwd")
      barlen <- nchar(orient.fwd)
      fq.rc <- reverseComplement(fq)
      keepF <- narrow(sread(fq), 1, barlen) == orient.fwd
      keepR <- narrow(sread(fq.rc), 1, barlen) == orient.fwd & 
        !keepF
      fq <- ShortReadQ(sread = c(sread(fq[keepF]), sread(fq.rc[keepR])), 
                       quality = c(quality(quality(fq[keepF])), quality(quality(fq.rc[keepR]))), 
                       id = c(id(fq[keepF]), id(fq.rc[keepR])))
    }
    if (is.finite(maxLen)) {
      fq <- fq[width(fq) <= maxLen]
    }
    fq <- fq[width(fq) >= start]
    fq <- narrow(fq, start = start, end = NA)
    if (trimRight > 0) {
      fq <- fq[width(fq) > trimRight]
      fq <- narrow(fq, start = NA, end = width(fq) - trimRight)
    }
    enc <- encoding(quality(fq))
    if (is.numeric(truncQ)) {
      ind <- which(enc == truncQ)
      if (length(ind) != 1) 
        stop("Encoding for this truncQ value not found.")
      truncQ <- names(enc)[[ind]]
    }
    if (length(fq) > 0) 
      fq <- trimTails(fq, 1, truncQ)
    truncQ <- enc[truncQ]
    if (!is.na(end)) {
      fq <- fq[width(fq) >= end]
    }
    fq <- narrow(fq, start = 1, end = end)
    fq <- fq[width(fq) >= minLen]
    fq <- fq[nFilter(maxN)(fq)]
    keep <- rep(TRUE, length(fq))
    qq <- as(quality(fq), "matrix")
    if (minQ > truncQ) 
      keep <- keep & (apply(qq, 1, min, na.rm = TRUE) > 
                        minQ)
    if (maxEE < Inf) {
      keep <- keep & C_matrixEE(qq) <= maxEE
    }
    fq <- fq[keep]
    if (rm.phix) {
      is.phi <- isPhiX(as(sread(fq), "character"), ...)
      fq <- fq[!is.phi]
    }
    if (rm.lowcomplex > 0) {
      sqcmplx <- seqComplexity(sread(fq), ...)
      fq <- fq[sqcmplx >= rm.lowcomplex]
    }
    outseqs <- outseqs + length(fq)
    if (first) {
      writeFastq(fq, fout, "w", compress = compress)
      first = FALSE
    } else {
      writeFastq(fq, fout, "a", compress = compress)
    }
  }
  if (verbose) {
    outperc <- round(outseqs * 100/inseqs, 1)
    outperc <- paste(" (", outperc, "%)", sep = "")
    message("Read in ", inseqs, ", output ", outseqs, outperc, 
            " filtered sequences.", sep = "")
  }
  if (outseqs == 0) {
    message(paste("The filter removed all reads:", fout ))#, 
                # "not written."))
      if(!fout %in% list.files()) system(paste("touch",fout))
    #file.remove(fout)
  }
  return(invisible(c(reads.in = inseqs, reads.out = outseqs)))
}


# File parsing
fastq <- snakemake@input[[1]]
filt <- snakemake@output[[1]]
sampleName <- gsub("preprocessing/","run.",gsub("/","",gsub(".fastq","",fastq)))
filtpath <- gsub("/.+","",filt)

#output table with number of reads pre-filtering

print("filtering")

fastqFilter0(fastq, 
            filt,
            truncLen=snakemake@config[['filtering']][['trunc_length']][['fwd']], 
            maxEE=snakemake@config[['filtering']][['max_EE']][['fwd']],
            maxLen=snakemake@config[['filtering']][['maxLen']][['fwd']],
            minLen=snakemake@config[['filtering']][['minLen']][['fwd']],
            minQ=snakemake@config[['filtering']][['minQ']][['fwd']],
            truncQ=snakemake@config[['filtering']][['trunc_qual']], 
            maxN=snakemake@config[['filtering']][['maxN']],
            rm.phix=as.logical(snakemake@config[['filtering']][['rm_phix']]),
            compress=TRUE, verbose=TRUE)

#output table with number of reads post-filtering
print("filtering done")
