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
if(!require(dada2)){
  BiocManager::install("GenomeInfoDbData",update=F,ask=F)
  require(dada2)
}
library(ShortRead)

fastqPairedFilter0 <- function (fn, fout, maxN = c(0, 0), truncQ = c(2, 2), truncLen = c(0,0),
                                 maxLen = c(Inf, Inf), minLen = c(20, 20), trimLeft = c(0, 0), 
                                 trimRight = c(0, 0), minQ = c(0, 0), maxEE = c(Inf, Inf), 
                                 rm.phix = c(TRUE, TRUE), rm.lowcomplex = c(0, 0), matchIDs = FALSE, 
                                 orient.fwd = NULL, id.sep = "\\s", id.field = NULL, n = 1e+06, 
                                 OMP = TRUE, qualityType = "Auto", compress = TRUE, verbose = FALSE, ...) {
  if (!OMP) {
    ompthreads <- .Call(ShortRead:::.set_omp_threads, 1L)
    on.exit(.Call(ShortRead:::.set_omp_threads, ompthreads))
  }
  if (!is.character(fn) || length(fn) != 2) stop("Two paired input file names required.")
  if (!is.character(fout) || length(fout) != 2) stop("Two paired output file names required.")
  if (any(duplicated(c(fn, fout)))) {
    stop("The output and input file names must be different.")
  }
  for (var in c("maxN", "truncQ", "truncLen", "maxLen", "minLen", 
    "trimLeft", "trimRight", "minQ", "maxEE", "rm.phix", 
    "rm.lowcomplex")) {
    if (length(get(var)) == 1) {
      assign(var, c(get(var), get(var)))
    }
    if (length(get(var)) != 2) {
      stop(paste("Input variable", var, "must be length 1 or 2 (Forward, Reverse)."))
    }
  }
  startF <- max(1, trimLeft[[1]] + 1, na.rm = TRUE)
  startR <- max(1, trimLeft[[2]] + 1, na.rm = TRUE)
  endF <- truncLen[[1]]
  if (endF < startF) {
    endF = NA
  }
  endF <- endF - startF + 1
  endR <- truncLen[[2]]
  if (endR < startR) {
    endR = NA
  }
  endR <- endR - startR + 1
  fF <- FastqStreamer(fn[[1]], n = n)
  on.exit(close(fF))
  fR <- FastqStreamer(fn[[2]], n = n)
  on.exit(close(fR), add = TRUE)
  if (file.exists(fout[[1]])) {
    if (file.remove(fout[[1]])) {
      if (verbose) message("Overwriting file:", fout[[1]])
    } else {
      stop("Failed to overwrite file:", fout[[1]])
    }
  }
  if (file.exists(fout[[2]])) {
    if (file.remove(fout[[2]])) {
      if (verbose) message("Overwriting file:", fout[[2]])
    } else {
      stop("Failed to overwrite file:", fout[[2]])
    }
  }
  first = TRUE
  remainderF <- ShortReadQ()
  remainderR <- ShortReadQ()
  casava <- "Undetermined"
  inseqs = 0
  outseqs = 0
  while (TRUE) {
    suppressWarnings(fqF <- yield(fF, qualityType = qualityType))
    suppressWarnings(fqR <- yield(fR, qualityType = qualityType))
    if (length(fqF) == 0 && length(fqR) == 0) {
      break
    }
    inseqs <- inseqs + length(fqF)
    if (matchIDs) {
      if (first) {
        if (is.null(id.field)) {
          id1 <- as.character(id(fqF)[[1]])
          id.fields <- strsplit(id1, id.sep)[[1]]
          ncolon <- sapply(gregexpr(":", id.fields), 
            length)
          ncoltab <- table(ncolon)
          if (max(ncolon) == 6 && ncoltab["6"] == 1) {
            casava <- "Current"
            id.field <- which(ncolon == 6)
          } else if (max(ncolon) == 4 && ncoltab["4"] == 
            1) {
            casava <- "Old"
            id.field <- which(ncolon == 4)
          } else {
            stop("Couldn't automatically detect the sequence identifier field in the fastq id string.")
          }
        }
      } else {
        fqF <- append(remainderF, fqF)
        fqR <- append(remainderR, fqR)
      }
    } else {
      if (length(fqF) != length(fqR)) 
        stop("Mismatched forward and reverse sequence files: ", 
          length(fqF), ", ", length(fqR), ".")
    }
    if (matchIDs) {
      idsF <- sapply(strsplit(as.character(id(fqF)), id.sep), 
        `[`, id.field)
      idsR <- sapply(strsplit(as.character(id(fqR)), id.sep), 
        `[`, id.field)
      if (casava == "Old") {
        idsF <- sapply(strsplit(idsF, "#"), `[`, 1)
      }
      lastF <- max(c(0, which(idsF %in% idsR)))
      lastR <- max(c(0, which(idsR %in% idsF)))
      if (lastF < length(fqF)) {
        remainderF <- fqF[(lastF + 1):length(fqF)]
      } else {
        remainderF <- ShortReadQ()
      }
      if (lastR < length(fqR)) {
        remainderR <- fqR[(lastR + 1):length(fqR)]
      } else {
        remainderR <- ShortReadQ()
      }
      fqF <- fqF[idsF %in% idsR]
      fqR <- fqR[idsR %in% idsF]
    }
    if (!is.null(orient.fwd)) {
      if (!C_isACGT(orient.fwd)) 
        stop("Non-ACGT characters detected in orient.fwd")
      barlen <- nchar(orient.fwd)
      keepF <- narrow(sread(fqF), 1, barlen) == orient.fwd
      keepR <- (narrow(sread(fqR), 1, barlen) == orient.fwd) & 
        !keepF
      fq <- ShortReadQ(sread = c(sread(fqF[keepF]), sread(fqR[keepR])), 
        quality = c(quality(quality(fqF[keepF])), quality(quality(fqR[keepR]))), 
        id = c(id(fqF[keepF]), id(fqR[keepR])))
      fqR <- ShortReadQ(sread = c(sread(fqR[keepF]), sread(fqF[keepR])), 
        quality = c(quality(quality(fqR[keepF])), quality(quality(fqF[keepR]))), 
        id = c(id(fqR[keepF]), id(fqF[keepR])))
      fqF <- fq
      rm(fq)
    }
    if (is.finite(maxLen[[1]]) || is.finite(maxLen[[2]])) {
      keep <- width(fqF) <= maxLen[[1]] & width(fqR) <= 
        maxLen[[2]]
      fqF <- fqF[keep]
      fqR <- fqR[keep]
    }
    keep <- (width(fqF) >= startF & width(fqR) >= startR)
    fqF <- fqF[keep]
    fqF <- narrow(fqF, start = startF, end = NA)
    fqR <- fqR[keep]
    fqR <- narrow(fqR, start = startR, end = NA)
    if (trimRight[[1]] > 0) {
      keep <- width(fqF) > trimRight[[1]]
      fqF <- fqF[keep]
      fqR <- fqR[keep]
      fqF <- narrow(fqF, start = NA, end = width(fqF) - 
        trimRight[[1]])
    }
    if (trimRight[[2]] > 0) {
      keep <- width(fqR) > trimRight[[2]]
      fqF <- fqF[keep]
      fqR <- fqR[keep]
      fqR <- narrow(fqR, start = NA, end = width(fqR) - 
        trimRight[[2]])
    }
    encF <- encoding(quality(fqF))
    encR <- encoding(quality(fqR))
    if (is.numeric(truncQ)) {
      indF <- which(encF == truncQ[[1]])
      indR <- which(encR == truncQ[[2]])
      if (!(length(indF) == 1 && length(indR) == 1)) stop("Encoding for this truncQ value not found.")
      truncQ <- c(names(encF)[[indF]], names(encR)[[indR]])
    }
    if (length(fqF) > 0) {
      rngF <- trimTails(fqF, 1, truncQ[[1]], ranges = TRUE)
      fqF <- narrow(fqF, 1, end(rngF))
    }
    if (length(fqR) > 0) {
      rngR <- trimTails(fqR, 1, truncQ[[2]], ranges = TRUE)
      fqR <- narrow(fqR, 1, end(rngR))
    }
    truncQ <- c(encF[truncQ[1]], encR[truncQ[2]])
    keep <- (width(fqF) > 0 & width(fqR) > 0)
    fqF <- fqF[keep]
    fqR <- fqR[keep]
    keep <- rep(TRUE, length(fqF))
    if (!is.na(endF)) {
      keep <- keep & (width(fqF) >= endF)
    }
    if (!is.na(endR)) {
      keep <- keep & (width(fqR) >= endR)
    }
    fqF <- fqF[keep]
    fqR <- fqR[keep]
    fqF <- narrow(fqF, start = 1, end = endF)
    fqR <- narrow(fqR, start = 1, end = endR)
    keep <- width(fqF) >= minLen[[1]] & width(fqR) >= minLen[[2]]
    fqF <- fqF[keep]
    fqR <- fqR[keep]
    suppressWarnings(keep <- nFilter(maxN[[1]])(fqF) & nFilter(maxN[[2]])(fqR))
    fqF <- fqF[keep]
    fqR <- fqR[keep]
    keep <- rep(TRUE, length(fqF))
    qmat <- as(quality(fqF), "matrix")
    if (minQ[[1]] > truncQ[[1]]) suppressWarnings(keep <- keep & (apply(qmat, 1, 
        min, na.rm = TRUE) > minQ[[1]]))
    if (maxEE[[1]] < Inf) keep <- keep &  .Call('_dada2_C_matrixEE', PACKAGE = 'dada2',qmat) <= maxEE[[1]]
    qmat <- as(quality(fqR), "matrix")
    if (minQ[[2]] > truncQ[[2]]) suppressWarnings(keep <- keep & (apply(qmat, 1, 
        min, na.rm = TRUE) > minQ[[2]]))
    if (maxEE[[2]] < Inf) keep <- keep &  .Call('_dada2_C_matrixEE', PACKAGE = 'dada2',qmat) <= maxEE[[2]]
    fqF <- fqF[keep]
    fqR <- fqR[keep]
    rm(qmat)
    if (length(fqF) != length(fqR)) stop("Filtering caused mismatch between forward and reverse sequence lists: ", 
        length(fqF), ", ", length(fqR), ".")
    if (rm.phix[[1]] && rm.phix[[2]]) {
      is.phi <- isPhiX(as(sread(fqF), "character"), ...)
      is.phi <- is.phi | isPhiX(as(sread(fqR), "character"), 
        ...)
    } else if (rm.phix[[1]] && !rm.phix[[2]]) {
      is.phi <- isPhiX(as(sread(fqF), "character"), ...)
    } else if (!rm.phix[[1]] && rm.phix[[2]]) {
      is.phi <- isPhiX(as(sread(fqR), "character"), ...)
    }
    if (any(rm.phix)) {
      fqF <- fqF[!is.phi]
      fqR <- fqR[!is.phi]
    }
    if (rm.lowcomplex[[1]] > 0 && rm.lowcomplex[[2]] > 0) {
      is.lowc <- (seqComplexity(sread(fqF), ...) < rm.lowcomplex[[1]])
      is.lowc <- is.lowc | (seqComplexity(sread(fqF), 
        ...) < rm.lowcomplex[[2]])
    } else if (rm.lowcomplex[[1]] && !rm.lowcomplex[[2]]) {
      is.lowc <- (seqComplexity(sread(fqF), ...) < rm.lowcomplex[[1]])
    } else if (!rm.lowcomplex[[1]] && rm.lowcomplex[[2]]) {
      is.lowc <- (seqComplexity(sread(fqR), ...) < rm.lowcomplex[[2]])
    }
    if (rm.lowcomplex[[1]] > 0 && rm.lowcomplex[[2]] > 0) {
      fqF <- fqF[!is.lowc]
      fqR <- fqR[!is.lowc]
    }
    outseqs <- outseqs + length(fqF)
    if (first) {
      writeFastq(fqF, fout[[1]], "w", compress = compress)
      writeFastq(fqR, fout[[2]], "w", compress = compress)
      first = FALSE
    }
    else {
      writeFastq(fqF, fout[[1]], "a", compress = compress)
      writeFastq(fqR, fout[[2]], "a", compress = compress)
    }
  }
#  if (outseqs == 0) {
#  }
  if (verbose) {
    outperc <- round(outseqs * 100/inseqs, 1)
    outperc <- paste(" (", outperc, "%)", sep = "")
    message("Read in ", inseqs, " paired-sequences, output ", 
      outseqs, outperc, " filtered paired-sequences.", 
      sep = "")
  }
  if (outseqs == 0) {
    message(paste("The filter removed all reads:", fout[[1]], 
      "and", fout[[2]] ))#, "not written."))
    if(!fout[[1]] %in% list.files()) system(paste("touch",fout[[1]])) 
    if(!fout[[2]] %in% list.files()) system(paste("touch",fout[[2]]))
 #   file.remove(fout[[1]])
 #   file.remove(fout[[2]])
  }
  return(invisible(c(reads.in = inseqs, reads.out = outseqs)))
}



# File parsing
fastqF <- snakemake@input[[1]]
fastqR <- snakemake@input[[2]]
filtF <- snakemake@output[[1]]
filtR <- snakemake@output[[2]]
sampleName <- gsub("preprocessing/","run.",gsub("/","",gsub(".fwd.fastq","",fastqF)))
#filtpath <- gsub("/.+","",filtF)

#output table with number of reads pre-filtering

print(paste("filtering",sampleName))

print(c(snakemake@config[['filtering']][['trunc_qual']][['fwd']],
                               snakemake@config[['filtering']][['trunc_qual']][['rvs']]))

fastqPairedFilter0(c(fastqF,fastqR), 
                    c(filtF,filtR),
                    truncLen=c(snakemake@config[['filtering']][['trunc_length']][['fwd']],
                               snakemake@config[['filtering']][['trunc_length']][['rvs']]), 
                    maxEE=c(snakemake@config[['filtering']][['max_EE']][['fwd']],
                               snakemake@config[['filtering']][['max_EE']][['rvs']]), 
                    maxLen=c(snakemake@config[['filtering']][['maxLen']][['fwd']],
                               snakemake@config[['filtering']][['maxLen']][['rvs']]),       
                    minLen=c(snakemake@config[['filtering']][['minLen']][['fwd']],
                               snakemake@config[['filtering']][['minLen']][['rvs']]),
                    minQ=c(snakemake@config[['filtering']][['minQ']][['fwd']],
                               snakemake@config[['filtering']][['minQ']][['rvs']]),
                    truncQ=c(snakemake@config[['filtering']][['trunc_qual']][['fwd']],
                               snakemake@config[['filtering']][['trunc_qual']][['rvs']]),
                    maxN=snakemake@config[['filtering']][['maxN']],
                    rm.phix=as.logical(snakemake@config[['filtering']][['rm_phix']]),
                    trimLeft=c(snakemake@config[['filtering']][['trim_left']][['fwd']],
                               snakemake@config[['filtering']][['trim_left']][['rvs']]), 
                    compress=TRUE, verbose=TRUE)

#output table with number of reads post-filtering
print("filtering done")
